# Documentation source open questions

Date: 2026-07-30

## Completed review gates

`@/workflows/review-documentation-sources.donna.md` reviewed these complete
sources:

- `@/R/ProBatchFeatures.R`
- `@/R/correct_batch_effects.R`
- `@/R/matrix_adapter.R`
- `@/R/pb_missing_filters.R`
- `@/R/plot_split_violinplot.R`
- `@/R/proteome_wide_diagnostics.R`
- `@/R/registry.R`
- `@/R/step_result.R`

All configured Depmesh relations were queried separately for every source.
The deterministic workflow gate confirmed that the sources parse, relation
output contains no generated files, and changes since the saved workflow
baseline affect only Roxygen2 lines. No command in this workflow read,
enumerated, compared, generated, or modified `man/` or `NAMESPACE`; no
documentation-generation, package-wide, network, staging, commit, or session
reset action occurred.

After the four source repairs below, the final eight-source review was
repeated. Its workflow artifact validation, rendering, scope validation, and
deterministic documentation-source gate all passed. The final gate again
confirmed R parsing, generated-file-free relation checks, and Roxygen-only
source-integrity enforcement. No Roxygen edit was needed in that run; semantic
completion remains blocked only by the t-SNE product decision recorded below.

## Source behavior repairs completed outside the documentation workflow

The documentation-source workflow cannot change package behavior. The semantic
review identified the following small implementation repairs. They were
completed after that workflow closed and their four focused test contexts
passed without unexpected warnings:

1. Source: `@/R/proteome_wide_diagnostics.R`
   Public interface: `plot_PCA()`
   Evidence: the shared embedding documentation requires a finite numeric
   `fill_the_missing` value, and t-SNE/UMAP call
   `.pb_validate_embedding_missing()`, but PCA does not. This makes PCA accept
   `Inf` even though the three APIs are required to share the same policy.
   Resolution: PCA now applies the shared validator, with focused diagnostic
   coverage.

2. Source: `@/R/matrix_adapter.R`
   Public interface: `pb_apply_matrix_method()`
   Evidence: the public contract and
   `@/specs/behavior/core_baseline.md` define only `error`, `keep`,
   `drop_features`, and `fill`, but the shared compatibility normalizer also
   accepts deprecated correction-only aliases.
   Resolution: the adapter now rejects non-canonical policies, with focused
   adapter coverage.

3. Source: `@/R/correct_batch_effects.R`
   Public interface: `adjust_batch_trend_dm()`
   Evidence: the function declares `...`, and the shared documentation states
   that trend-adjustment arguments reach the selected fitter, but the matrix
   wrapper currently drops them when it calls `adjust_batch_trend_df()`.
   Resolution: the matrix wrapper now forwards `...`, with a small direct
   forwarding test.

4. Source: `@/R/pb_missing_filters.R`
   Public interface: `pb_filterNA()` and `pb_groupfilterNA()`
   Evidence: a scalar explicit `final_name` is correctly rejected for
   multi-assay input, but the error still says that length one is accepted.
   Resolution: the error now states the enforced one-name-per-assay rule, and
   existing focused coverage asserts it.

## Maintainer decision

Source: `@/R/proteome_wide_diagnostics.R`

Public interface: `plot_TSNE()`

Evidence: the baseline requires at least two samples and clamps perplexity to
`max(1, floor((n_samples - 1) / 3))`, while the pinned optional backend
`Rtsne 0.17` rejects two- and three-sample inputs at perplexity one. The
reproduction and alternatives are also recorded in
`@/workflows/reports/open-questions.md`.

Ownership question: whether Core should require four samples for the Rtsne
backend or specify a different small-sample behavior.

Required decision: choose and specify the supported small-sample contract
before changing source documentation, behavior, or tests. After all source
repairs, the maintainer must manually regenerate documentation and `NAMESPACE`
and inspect the installed public surface.

No cataloged workflow owns this product decision. Once the contract is chosen
and specified, the affected implementation and test work should use the normal
specification-verification and focused-test workflows before this
documentation-source review is rerun.
