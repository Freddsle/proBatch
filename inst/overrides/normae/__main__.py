#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# NormAE CLI overlay for proBatch integration (without PCA plotting)
# From https://github.com/luyiyun/NormAE

import argparse
import os
import sys
import pandas as pd
import numpy as np
import importlib

def _none_if_empty(s):
    if s is None:
        return None
    s = str(s).strip()
    if s in ("", "''", '""', "None", "none", "NA", "na"):
        return None
    return s

def _coerce_float_matrix(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if out.isnull().values.any():
        # NormAE expects no NaNs; caller (proBatch) handles imputation/removal already.
        raise SystemExit("ERROR: NaNs detected in meta_csv after numeric coercion.")
    return out

def _infer_id_col(samples: pd.DataFrame, X_idx: pd.Index) -> str:
    # pick the column with maximal exact-match coverage to X index
    Xs = set(map(str, X_idx.tolist()))
    best = None
    best_hits = -1
    for col in samples.columns:
        vals = list(map(str, samples[col].tolist()))
        hits = sum(v in Xs for v in vals)
        if hits > best_hits and len(set(vals)) == len(vals):  # uniqueness safeguard
            best = col
            best_hits = hits
    if best is None:
        # fallback to first column
        return samples.columns[0]
    return best

def _mask_equals(values: pd.Series, target: str) -> pd.Series:
    # robust QC match
    v = values.astype(str).str.strip().str.lower()
    t = str(target).strip().lower()
    return v == t

def main():
    parser = argparse.ArgumentParser(
        description="NormAE CLI overlay without PCA plotting"
    )
    parser.add_argument("--meta_csv", required=True, help="CSV with features (columns) and samples (rows).")
    parser.add_argument("--sample_csv", required=True, help="CSV with sample annotation.")
    parser.add_argument("--output_dir", required=True, help="Output directory containing X_clean.csv")
    parser.add_argument("--batch_indicator_col", required=True, help="Column name in sample_csv for batch labels.")
    parser.add_argument("--order_indicator_col", default="", help="Column name in sample_csv for injection order.")
    parser.add_argument("--qc_indicator_col", default="", help="Column in sample_csv marking QC samples.")
    parser.add_argument("--qc_indicator_value", default="QC", help="Value indicating QC rows when qc_indicator_col is provided.")

    # Accept and ignore any additional CLI args (keeps compatibility with future flags).
    args, _extras = parser.parse_known_args()

    order_col = _none_if_empty(args.order_indicator_col)
    qc_col    = _none_if_empty(args.qc_indicator_col)
    qc_val    = args.qc_indicator_value

    # 1) Load data
    X = pd.read_csv(args.meta_csv, index_col=0)
    # proBatch writes meta as t(data_matrix): rows=samples, cols=features
    X = _coerce_float_matrix(X)

    samples = pd.read_csv(args.sample_csv)
    if args.batch_indicator_col not in samples.columns:
        raise SystemExit(f"ERROR: batch column '{args.batch_indicator_col}' not found in sample_csv.")

    if order_col is not None and order_col not in samples.columns:
        raise SystemExit(f"ERROR: order column '{order_col}' not found in sample_csv.")

    id_col = _infer_id_col(samples, X.index)
    samples = samples.set_index(id_col)

    # Strict reindex to X order; error if missing
    try:
        samples = samples.loc[X.index]
    except KeyError as e:
        missing = [s for s in X.index if s not in samples.index]
        raise SystemExit(f"ERROR: sample_csv is missing {len(missing)} samples present in meta_csv; first missing: {missing[:5]}")

    # 2) Prepare labels / orders / QC
    y = samples[args.batch_indicator_col].values
    z = samples[order_col].values if order_col is not None else None

    X_qc = y_qc = z_qc = None
    if qc_col is not None:
        if qc_col not in samples.columns:
            raise SystemExit(f"ERROR: QC column '{qc_col}' not found in sample_csv.")
        qc_mask = _mask_equals(samples[qc_col], qc_val)
        if qc_mask.any():
            X_qc = X.loc[qc_mask].values
            y_qc = samples.loc[qc_mask, args.batch_indicator_col].values
            z_qc = samples.loc[qc_mask, order_col].values if order_col is not None else None

    # 3) Import the real NormAE class from the installed package
    #    (we are only overriding __main__, not the rest of the package).
    try:
        normae_pkg = importlib.import_module("normae")
        if hasattr(normae_pkg, "NormAE"):
            NormAE = normae_pkg.NormAE
        else:
            # Some layouts expose class under submodule:
            NormAE = importlib.import_module("normae.normae").NormAE
    except Exception as e:
        raise SystemExit(f"ERROR: failed to import NormAE from installed package: {e}")

    # 4) Fit / Transform (no plotting)
    model = NormAE()
    model.fit(X.values, y=y, z=z, X_qc=X_qc, y_qc=y_qc, z_qc=z_qc)
    X_clean = model.transform(X.values)

    # 5) Write output exactly where proBatch expects it
    os.makedirs(args.output_dir, exist_ok=True)
    out_csv = os.path.join(args.output_dir, "X_clean.csv")
    pd.DataFrame(X_clean, index=X.index, columns=X.columns).to_csv(out_csv)

if __name__ == "__main__":
    main()
