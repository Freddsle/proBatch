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

## 2026-08-01 clean-export review

The clean-export documentation review validated and rendered the workflow,
captured a protected baseline for 26 R sources, read every scoped source, and
queried all configured Depmesh relations separately for each source. The
review repaired only Roxygen2 lines, including stale packaged-data dimensions
and classes, argument/default descriptions, conditional and multi-assay return
values, and examples. It did not inspect or modify `man/` or `NAMESPACE` and
did not generate documentation, run examples, stage files, create commits, or
change package behavior.

The deterministic post-edit gate was not run because semantic inspection found
the behavior-level questions below and the workflow therefore followed its
required blocker path.

### Ignored base-size argument

Source: `@/R/feature_level_diagnostics.R`

Public interface: `plot_peptides_of_one_protein()`

Evidence: the public function exposes and inherits documentation for
`base_size`, but its call to `plot_single_feature()` omits that argument. A
non-default value is accepted and silently ignored.

Ownership question: whether the function should honor its existing public
argument or remove/deprecate that argument and document the compatibility
path.

Required decision: authorize the intended public contract. If `base_size`
remains supported, forward it and add focused regression coverage.

### Custom measurement column omitted during combined normalization

Source: `@/R/normalize.R`

Public interface: `normalize_data_df()`

Evidence: the public function accepts `measure_col`, but when `log_base` is not
`NULL` its call to `log_transform_df()` omits `measure_col`. The optional log
step therefore operates on the default `Intensity` column instead of the
requested measurement column.

Ownership question: whether the combined wrapper is required to honor the
existing custom-column contract during every enabled step.

Required decision: authorize forwarding `measure_col` to the log step and add
focused coverage, or explicitly narrow and deprecate the current public
argument contract.

### Custom feature identifier omitted from median normalization

Source: `@/R/normalize.R`

Public interface: `normalize_data_df()`

Evidence: the median-centering branch calls `normalize_sample_medians_df()`
without forwarding the public `feature_id_col` argument. Non-default feature
identifier columns therefore do not receive the contract exposed by the
wrapper.

Ownership question: whether the wrapper should consistently propagate custom
identifier columns to its selected implementation.

Required decision: authorize forwarding `feature_id_col` and add focused
regression coverage, or explicitly narrow and deprecate the wrapper contract.

No cataloged workflow can make these product decisions. If behavior repairs
are authorized, their specification impact should be assessed with
`@/workflows/verify-specifications.donna.md`, their tests run with
`@/workflows/run-tests.donna.md`, and this documentation-source workflow rerun
before the maintainer regenerates documentation and `NAMESPACE`.
