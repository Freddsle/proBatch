# Open questions

## Minimum sample count for t-SNE

Status: open

Raised: 2026-07-30

`./specs/behavior/core_baseline.md` requires `plot_TSNE()` to accept at least
two samples and to reduce excessive perplexity to
`max(1, floor((n_samples - 1) / 3))`.

That contract is not executable with the declared backend. With `Rtsne` 0.17,
the required perplexity is `1` for both two and three samples, and
`Rtsne::Rtsne()` rejects both inputs with
`perplexity is too large for the number of samples`. Four samples succeed with
the same perplexity.

A maintainer decision is required before changing public behavior:

- raise the minimum supported sample count to four and align the specification,
  implementation, documentation, and tests; or
- specify a different two- and three-sample behavior that does not claim to use
  `Rtsne` for those inputs.

Until that decision is made, the implementation and specification remain
unchanged.

## Complete-suite thread creation in the Donna Pixi environment

Status: resolved by maintainer verification

Raised: 2026-07-30

Reconfirmed: 2026-07-31

Resolved: 2026-07-31

In the 2026-07-30 run, the final Donna repair boundary contained
`@/R/plot_helpers.R` and
`@/tests/testthat/test-proteome_wide_diagnostics.R`. Workflow artifact,
test-environment, and plan validation passed, as did the selected context:

- `test-proteome_wide_diagnostics.R`

The immediately preceding focused run passed `test-plot_helpers.R` and
`test-proteome_wide_diagnostics.R`. The preceding compatibility run also passed
`test-ProBatchFeatures.R`, `test-batch_effect_steps.R`,
`test-design_diagnostics.R`, and `test-matrix_adapter.R`. Earlier focused runs
of the correction and missing-data contexts passed as recorded in Donna's
completed history.

The final complete source-level testthat run passed every context except the
pre-existing normalization context, which reported these two errors:

```text
test-normalize.R:5:5
Error in normalize.quantiles(data_matrix):
ERROR; return code from pthread_create() is 22

test-normalize.R:110:5
Error in normalize.quantiles(data_matrix):
ERROR; return code from pthread_create() is 22
```

Both failures originate in `preprocessCore::normalize.quantiles()` through
`@/R/normalize.R`. That source and `@/tests/testthat/test-normalize.R` are
outside the selected repair boundary and were unchanged. The process also
reported restricted system-bus/`timedatectl` diagnostics, which did not become
test failures.

The maintainer reports that `devtools::test()` succeeds without warnings or
errors in the normal development environment.

The final 2026-07-31 check-cleanup test workflow reconfirmed the same
limitation. Its authorized boundary contained nine R sources, five changed
focused tests, and `inst/WORDLIST`. Donna's eleven selected contexts all
passed, including the portable plot-label rendering, PCA dots contract,
data.table selection, plotting-helper, and batch-correction cases. The
required complete run then passed every other context and reproduced only the
same two `test-normalize.R` errors shown above.

An earlier reproduction with sandbox escalation produced the identical
`preprocessCore::normalize.quantiles()` error. This rules out the cleanup diff
and ordinary sandbox permissions; the precise cause within the existing
Pixi/conda threading runtime remains undiagnosed. Each complete attempt removed
the newly generated `tests/testthat/Rplots.pdf` and preserved the governed
worktree, `pixi.toml`, and `pixi.lock`.

No remaining cataloged workflow is authorized to provision or rebuild
dependencies, and the test workflow explicitly forbids doing so. The
maintainer's subsequent `devtools::check()` completed the package tests and
did not reproduce either normalization failure, confirming that this blocker
was confined to the agent's Pixi environment. No normalization or test change
is justified by the environment-only failure. That rerun exposed only a
separate spelling NOTE, whose two Roxygen source phrases were then corrected.
