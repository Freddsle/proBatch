#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NormAE CLI overlay for proBatch integration (no plotting).

Goal
----
Behave EXACTLY like upstream `python -m normae` while preventing any plotting
and avoiding GUI/file I/O from matplotlib, without re-implementing CLI logic.

Strategy
--------
1) Use non-GUI backend; no-op `plt.show`, `plt.savefig`, and Figure.savefig.
2) Patch NormAE's own plotters BEFORE the CLI binds them:
     - `normae.utils.plot_pca`  -> returns (Figure, Axes)
     - `normae.estimator.NormAE.plot_history` -> no-op (returns None)
3) Execute the official package entrypoint via `runpy.run_module("normae", "__main__")`.

Non-zero SystemExit codes propagate to the caller (e.g., your R wrapper).
"""

from __future__ import annotations

import os
import sys
import importlib
import runpy
from types import ModuleType
from typing import Callable

# --------------------------- generic helpers ---------------------------

def _noop(*args, **kwargs):
    return None

def _ensure_headless_matplotlib() -> None:
    """Select a headless backend and suppress display/writes.
    DO NOT patch plt.subplots/plt.figure — upstream relies on their return values.
    """
    os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        import matplotlib
        try:
            matplotlib.use("Agg", force=True)
        except Exception:
            pass

        # Patch pyplot functions that cause UI or disk writes
        try:
            import matplotlib.pyplot as plt
            if hasattr(plt, "show"):
                plt.show = _noop
            if hasattr(plt, "savefig"):
                plt.savefig = _noop
        except Exception:
            pass

        # Patch Figure.savefig too (covers fig.savefig(...))
        try:
            from matplotlib.figure import Figure
            if hasattr(Figure, "savefig"):
                Figure.savefig = _noop
            # Optional belt-and-suspenders: avoid layout churn
            if hasattr(Figure, "tight_layout"):
                # harmless if left alone, but neutralize just in case
                Figure.tight_layout = (lambda self, *a, **k: None)
        except Exception:
            pass

    except Exception:
        # If matplotlib is not installed, it's fine because we also patch NormAE plotters.
        pass


# --------------------------- NormAE patching ---------------------------

def _patch_callables_with_name_fragment(mod: ModuleType, frag: str, replacer: Callable) -> int:
    """Replace any callable attribute whose name contains `frag` (case-insensitive)."""
    n = 0
    frag = frag.lower()
    for attr in list(dir(mod)):
        if frag in attr.lower():
            try:
                obj = getattr(mod, attr)
            except Exception:
                continue
            if callable(obj):
                try:
                    setattr(mod, attr, replacer)
                    n += 1
                except Exception:
                    pass
    return n

def _plot_pca_stub(*args, **kwargs):
    """Return a real (Figure, Axes) pair so upstream `fig, _ = plot_pca(...)` keeps working."""
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1)  # returns (Figure, Axes) by contract
        return fig, ax
    except Exception:
        # Fallback dummy with the minimal methods sometimes used by callers
        class _Dummy:
            def savefig(self, *a, **k): return None
            def tight_layout(self, *a, **k): return None
        return _Dummy(), _Dummy()

def _plot_history_stub(self, *args, **kwargs):
    """Harmless replacement for NormAE.plot_history; upstream doesn't rely on its return."""
    return None

def _patch_normae_plotters() -> None:
    """Import NormAE components and neuter plotting *before* running the CLI."""
    # Import package root so submodules resolve
    importlib.import_module("normae")

    # 1) Patch utils.plot_pca specifically (must return (fig, ax))
    try:
        utils = importlib.import_module("normae.utils")
        if hasattr(utils, "plot_pca"):
            setattr(utils, "plot_pca", _plot_pca_stub)
        # For any other utils.*plot* helpers that might be called, default to true no-op.
        # (If the CLI later unpacks another function, we can specialize it similarly.)
        _patch_callables_with_name_fragment(utils, "plot", _noop)
        # Reassert the specific plot_pca stub in case the blanket patch hit it:
        if hasattr(utils, "plot_pca"):
            setattr(utils, "plot_pca", _plot_pca_stub)
    except Exception:
        pass

    # 2) Patch estimator.NormAE.plot_history
    try:
        estimator = importlib.import_module("normae.estimator")
        NormAE_cls = getattr(estimator, "NormAE", None)
        if NormAE_cls is not None and hasattr(NormAE_cls, "plot_history"):
            setattr(NormAE_cls, "plot_history", _plot_history_stub)
    except Exception:
        pass


# ------------------------------ main ------------------------------

def main() -> None:
    _ensure_headless_matplotlib()
    _patch_normae_plotters()

    # Execute upstream CLI exactly like `python -m normae`
    try:
        if sys.argv:
            sys.argv[0] = "normae"
        runpy.run_module("normae", run_name="__main__", alter_sys=True)
    except SystemExit as e:
        # propagate non-zero exits
        code = e.code
        if code not in (None, 0):
            raise

if __name__ == "__main__":
    main()
