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

## Complete-suite thread creation in the agent sandbox

Status: environment-limited

Raised: 2026-07-30

The final Donna repair boundary contained `@/R/plot_helpers.R` and
`@/tests/testthat/test-proteome_wide_diagnostics.R`. Workflow artifact
validation, test-environment validation, plan validation, and focused
execution passed. The selected context was:

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
errors in the normal development environment. Formal verification of the
exact post-audit tree therefore requires rerunning the complete suite in an
environment that can create the backend thread. No normalization or test
change is justified from this sandbox-only failure.
