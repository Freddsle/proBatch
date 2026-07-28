# BEC E. coli Core Migration Manifest

<!-- donna-migration scope=post-sync last-prefix=032 entries=31 workflows=24 commits=29 references=7 -->

## Immutable scope and ordering

- Source repository: `/home/yuliya/repos/other/proBatch`
- Source baseline: `ba6ee246eace090e71baa7aba302ca64e76ddb32`
- Synchronizing merge and foundation snapshot: `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`
- Source tip: `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`
- Stopped split repository: `/home/yuliya/repos/other/proBatch-core-split`
- Stopped split preflight: `49cee7cc978fbb149c262a5a783face32dd1d135`
- Stopped split implementation comparator: `29a7478dc7deea846a2c1ff1abd25a881e6f87db`
- Stopped split provenance note: `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md`
- Selected history: one foundation snapshot at the synchronizing merge, followed by the exact 29-commit output of `git rev-list --reverse --topo-order 5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab..e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.
- Ordering result: all 29 selected commits are single-parent and linear. The first parent is the foundation merge, each later parent is the preceding manifest commit, and the sequence ends at the pinned source tip.
- Exclusions: generated `man/**` and generated `NAMESPACE` are outside this audit and every child workflow. Source documentation means Roxygen2 comments under `R/`.

The synchronizing merge has parents `5b8579bc087396c94de2ed0da21fba34eea81106` and the named baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32`. Its date is `2026-03-12T16:12:32+01:00` and its complete message is `Merge remote-tracking branch 'origin/main' into BEC_ecoli_data`.

## Ownership boundary

The ownership map assigns 107 pre-split exports to core. The foundation contains 104 of them. Three core exports arrive later in the selected sequence: `plot_grouped_NA_heatmap` in commit 4, `pb_subset_samples` in commit 9, and `plot_NA_intensity` in commit 13.

The 26 companion-owned exports are:

- Provider correction: `correct_with_BERT`, `correct_with_NormAE`, `correct_with_PLSDA_batch`, `correct_with_RUVIII_C`, `correct_with_mComBat`, `correct_with_omicsGMF`, and `estimate_omicsGMF_rank`.
- Provider imputation: `imputeMissForest`, `imputeMissForest.ProBatchFeatures`, `imputeMissForest_df`, `imputeMissForest_dm`, `missForestImpute`, `imputePRONE`, `imputePRONE_df`, `imputePRONE_dm`, `impute_and_correct_with_omicsGMF`, `impute_with_omicsGMF`, and `impute_with_omicsGMF.ProBatchFeatures`.
- Benchmark metrics, embeddings, and plots: `calculate_classification_metrics`, `calculate_variance_partition`, `prepare_variance_partition_df`, `plot_variance_partition`, `plot_variance_partition.df`, `plot_intragroup_variation`, `plot_TSNE`, and `plot_UMAP`.

Provider implementations, provider registrations, benchmark orchestration, optional-provider discovery, Python integration, and companion-only dependencies remain outside core even when they share a source file or commit with a core hunk.

## Effective correction behavior at the source snapshot

The source package has no `Collate` field. In default C-locale source collation order, `R/correct_batch_effects_old.R` is sourced after `R/correct_batch_effects.R` and overwrites three differing definitions:

- `correct_batch_effects_df`
- `correct_batch_effects_dm`
- `correct_with_removeBatchEffect_dm`

The foundation and later children must compare the effective `_old.R` definitions for those symbols, not infer authority from a newer-looking filename. The stopped split consolidates correction definitions in `R/correct_batch_effects.R`; it is evidence for resolving the ambiguity, not an accepted patch.

## Stopped split evidence inventory

The stopped note says the implementation was rejected and must be reassessed from the baseline. The preflight changes only the `DESCRIPTION` mode from `100755` to `100644`. From the baseline through the implementation commit, excluding generated artifacts, the split changes 47 paths: 19 added and 28 modified.

- Project and package metadata: `.Rbuildignore`, `.gitignore`, `DESCRIPTION`, and `NEWS`.
- Modified R sources: `R/CV_calculation.R`, `R/ProBatchFeatures.R`, `R/auxiliary.R`, `R/correct_batch_effects.R`, `R/correlation-based_diagnostics.R`, `R/handle_missing_values.R`, `R/pb_missing_filters.R`, `R/plot_helpers.R`, `R/plot_missing.R`, `R/plot_split_violinplot.R`, `R/proBatch.R`, `R/proteome_wide_diagnostics.R`, `R/utility_funcs.R`, and `R/zzz_helpers.R`.
- Added R sources: `R/design_diagnostics.R`, `R/identifiers.R`, `R/matrix_adapter.R`, `R/metadata_diagnostics.R`, `R/pbf_input_helpers.R`, `R/plot_NA_intensity.R`, `R/registry.R`, and `R/step_result.R`.
- Installed resource: modified `inst/CITATION`.
- Added test helpers: `tests/testthat/helper-example-data.R` and `tests/testthat/helper-source-root.R`.
- Modified focused tests: `tests/testthat/test-ProBatchFeatures.R`, `test-auxiliary.R`, `test-batch_effect_steps.R`, `test-correct_batch_effect.R`, `test-correlation_based_diagnostics.R`, `test-handle_missing_values.R`, `test-pb_missing_helpers.R`, `test-proteome_wide_diagnostics.R`, and `test-utility_funcs.R`.
- Added focused tests: `tests/testthat/test-design_diagnostics.R`, `test-identifiers.R`, `test-lineage.R`, `test-matrix_adapter.R`, `test-metadata_diagnostics.R`, `test-plot_NA_intensity.R`, `test-registry.R`, `test-step_result.R`, and `test-symbol-ownership.R`.
- No editable data or vignette path changes in this baseline-to-implementation delta.

The split adds provider-aware registry, identifier, matrix-adapter, structured-result, and assay-lineage layers; consolidates duplicate definitions; removes companion providers; and hardens missing-value behavior. Every use below treats those changes as comparison evidence only.

## Ordered migration entries

### 002 — Foundation core API

<!-- donna-entry kind=foundation id=foundation-core-api slot=002 mode=workflow workflow=@/workflows/002_foundation_core_api.donna.md sha=- -->

- Snapshot: `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`, compared with baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32`.
- Classification: mixed foundation. Include all core-owned behavior at hunk level; exclude every companion provider, benchmark export, provider registration, benchmark script, optional-provider dependency, and companion-only vignette.
- Inclusion reason: the merge synchronizes the baseline into the accumulated BEC branch. A foundation is required to migrate the core API that predates the 29 post-synchronization commits without pretending the 282 older divergent commits were authored after synchronization.
- Core container and pipeline API from `R/ProBatchFeatures.R`: `ProBatchFeatures`, `ProBatchFeatures_from_long`, `as_ProBatchFeatures`, `pb_add_level`, `pb_aggregate_level`, `pb_as_long`, `pb_as_wide`, `pb_assay_matrix`, `pb_current_assay`, `pb_eval`, `pb_pipeline_name`, `get_chain`, `get_operation_log`, `pb_register_step`, `pb_list_steps`, `pb_has_step`, and `pb_transform`.
- Core conversion API from `R/auxiliary.R`: `long_to_matrix`, `matrix_to_long`, and `create_peptide_annotation`.
- Core missing-value API from `R/handle_missing_values.R` and `R/pb_missing_filters.R`: `handle_missing_values`, `pb_filterNA`, `pb_groupfilterNA`, `pb_infIsNA`, `pb_nNA`, and `pb_zeroIsNA`.
- Core metadata/design API from `R/utility_funcs.R`, `R/colors_for_annotation.R`, `R/date_conversion.R`, `R/design_diagnostics.R`, and `R/metadata_diagnostics.R`: `check_sample_consistency`, `convert_annotation_classes`, `date_to_sample_order`, `dates_to_posix`, `define_sample_order`, `detect_nested_batches`, `detect_outlier_samples`, `filter_metadata_columns`, `find_duplicated_columns`, `guess_factor_columns_if_needed`, `handle_factor_numeric_overlap`, `merge_rare_levels`, `metadata_column_summary`, `subbatch_detection`, `summarize_design`, `validate_batch_design`, and `warn_unmapped_columns`.
- Core transformation API from `R/correct_batch_effects.R`, `R/correct_batch_effects_old.R`, `R/fit_non_linear.R`, `R/normalize.R`, and `R/transform_raw.R`: `adjust_batch_trend_df`, `adjust_batch_trend_dm`, `center_feature_batch`, `center_feature_batch_means_df`, `center_feature_batch_means_dm`, `center_feature_batch_medians_df`, `center_feature_batch_medians_dm`, `fit_nonlinear`, `log_transform_df`, `log_transform_dm`, `normalize_data_df`, `normalize_data_dm`, `normalize_sample_medians_df`, `normalize_sample_medians_dm`, `quantile_normalize_df`, `quantile_normalize_dm`, `unlog_df`, and `unlog_dm`.
- Core correction API from the effective correction sources: `correct_batch_effects`, `correct_batch_effects_df`, `correct_batch_effects_dm`, `correct_with_ComBat`, `correct_with_ComBat_df`, `correct_with_ComBat_dm`, `correct_with_removeBatchEffect`, `correct_with_removeBatchEffect_df`, and `correct_with_removeBatchEffect_dm`.
- Core numerical/PVCA API from `R/CV_calculation.R`, `R/correlation-based_diagnostics.R`, and the core hunks of `R/explained_variance_plots.R`: `calculate_feature_CV`, `calculate_peptide_corr_distr`, `calculate_sample_corr_distr`, `calculate_PVCA`, and `prepare_PVCA_df`.
- Core visualization API from `R/CV_calculation.R`, `R/plot_missing.R`, `R/proteome_wide_diagnostics.R`, `R/initial_assessment.R`, `R/correlation-based_diagnostics.R`, `R/feature_level_diagnostics.R`, `R/plot_split_violinplot.R`, and `R/colors_for_annotation.R`: `plot_CV_distr`, `plot_CV_distr.df`, `plot_NA_density`, `plot_NA_frequency`, `plot_NA_heatmap`, `plot_PCA`, `plot_PVCA`, `plot_PVCA.df`, `plot_PVCA_stacked_from_saved`, `plot_boxplot`, `plot_corr_matrix`, `plot_heatmap_diagnostic`, `plot_heatmap_generic`, `plot_hierarchical_clustering`, `plot_iRT`, `plot_peptide_corr_distribution`, `plot_peptide_corr_distribution.corrDF`, `plot_peptides_of_one_protein`, `plot_protein_corrplot`, `plot_sample_corr_distribution`, `plot_sample_corr_distribution.corrDF`, `plot_sample_corr_heatmap`, `plot_sample_mean`, `plot_single_feature`, `plot_spike_in`, `plot_split_violin_with_boxplot`, `plot_with_fitting_curve`, `generate_colors_for_numeric`, and `sample_annotation_to_colors`.
- Editable source/support targets: the core files named above plus `R/batch_correction_helpers.R`, core-only helpers in `R/pbf_input_helpers.R`, `R/plot_helpers.R`, `R/proBatch.R`, `R/utility_funcs.R`, and `R/zzz_helpers.R`. Preserve Roxygen comments with their core functions.
- Mixed-file exclusions: provider registration and PRONE behavior in `R/ProBatchFeatures.R`/`R/zzz_helpers.R`; mComBat branches in `R/correct_batch_effects.R`; variance-partition exports in `R/explained_variance_plots.R`; TSNE/UMAP in `R/proteome_wide_diagnostics.R`; provider/embedding/intragroup helpers in `R/plot_helpers.R`; and omicsGMF-only helpers in `R/utility_funcs.R`.
- Whole-file exclusions: `R/M_ComBat.R`, `R/calculate_intragroup.R`, `R/classification_metrics.R`, `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_NormAE.R`, `R/correct_batch_RUVIIIC.R`, `R/correct_batch_omicsGMF.R`, `R/impute_PRONE.R`, `R/impute_missForest.R`, `R/impute_omicsGMF.R`, `inst/overrides/**`, `inst/scripts/**`, companion tests, and `vignettes/batch_evaluation_metrics.Rmd`.
- Focused core tests: `test-CV_calculation.R`, `test-ProBatchFeatures.R`, `test-ProBatchFeatures_links.R`, `test-auxiliary.R`, core cases in `test-batch_effect_steps.R`, `test-colors_for_annotation.R`, `test-correct_batch_effect.R`, `test-correlation_based_diagnostics.R`, `test-date_conversion.R`, `test-design_diagnostics.R`, core PVCA cases in `test-explained_variance_plots.R`, `test-feature_level_diagnostics.R`, `test-fit_non_linear.R`, `test-handle_missing_values.R`, `test-initial_assessment.R`, `test-metadata_diagnostics.R`, `test-normalize.R`, `test-pb_missing_helpers.R`, core cases in `test-plot_helpers.R`, `test-plot_missing.R`, `test-plot_split_violinplot.R`, core cases in `test-proteome_wide_diagnostics.R`, `test-transform_raw.R`, and `test-utility_funcs.R`.
- Metadata/dependencies: derive the minimum core dependency set from actual retained code. Exclude PRONE, Rtsne, umap, variancePartition, BERT, PLSDAbatch, reticulate, omicsGMF, sgdGMF, SingleCellExperiment, and RUVIIIC unless independently required by a retained core hunk. Compare the split `DESCRIPTION` dependency cleanup, but do not copy its R-version or dependency choices without current-project review.
- Split comparison targets: `R/ProBatchFeatures.R`, `R/registry.R`, `R/identifiers.R`, `R/matrix_adapter.R`, `R/step_result.R`, `R/correct_batch_effects.R`, `R/handle_missing_values.R`, `R/design_diagnostics.R`, `R/metadata_diagnostics.R`, `R/pbf_input_helpers.R`, `R/plot_missing.R`, `R/plot_helpers.R`, `R/proteome_wide_diagnostics.R`, their focused tests, `DESCRIPTION`, and `NEWS`.

### 003 — PVCA implementation consolidation

<!-- donna-entry kind=commit id=post-sync-001-42544a21-pvca-impl-dedup slot=003 mode=workflow workflow=@/workflows/003_42544a21f10c_pvca_impl_dedup.donna.md sha=42544a21f10ca6960d3e4c44d2833f764054d721 -->

- Parent/date/message: parent `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`; `2026-03-13T18:12:02+01:00`; subject `Refactor PVCA functions to use implementation helpers for improved readability and maintainability`; body empty.
- Non-generated paths: `R/explained_variance_plots.R`, `R/proteome_wide_diagnostics.R`.
- Classification/disposition: core, included for semantic comparison.
- Hunks and core functions: captures six authoritative PVCA method bodies as `.pb_*_impl` aliases, replaces duplicate bodies with delegates, and makes `plot_PVCA.df` an S3 generic. Affects `calculate_PVCA.default`, `calculate_PVCA.ProBatchFeatures`, `plot_PVCA.default`, `plot_PVCA.ProBatchFeatures`, `prepare_PVCA_df.default`, `prepare_PVCA_df.ProBatchFeatures`, and `plot_PVCA.df`.
- Focused tests/dependencies: no test changed; verify PVCA and symbol-ownership coverage. Existing dependencies are `pvca`, `Biobase`, `ggplot2`, `dplyr`, and `gridExtra`.
- Split targets: `R/proteome_wide_diagnostics.R` and `tests/testthat/test-symbol-ownership.R`; the split consolidates the family and supersedes the alias/delegation mechanism.

### 004 — PRONE normalization integration

<!-- donna-entry kind=commit id=post-sync-002-db128668-prone-normalization-integration slot=004 mode=workflow workflow=@/workflows/004_db128668458a_prone_normalization_integration.donna.md sha=db128668458a58ad31f66be5b6e39e2fedadbbe1 -->

- Parent/date/subject: parent `42544a21f10ca6960d3e4c44d2833f764054d721`; `2026-03-14T21:06:55+01:00`; subject `Add PRONE normalization methods integration with proBatch`.
- Complete body:
  - `- Implemented .pb_register_prone_normalization_steps to register PRONE normalization methods with optional prefixes.`
  - `-  .pb_prone_normalization_method_names to retrieve and clean method names from PRONE.`
  - `- .pb_make_prone_norm_step_fun to generate normalization functions for each method.`
  - `- Added tests.`
  - `- Added a vignette .`
  - `- Updated zzz_helpers to automatically register PRONE normalization steps upon package load.`
- Non-generated paths: `DESCRIPTION`, `R/ProBatchFeatures.R`, `R/calculate_intragroup.R`, `R/impute_PRONE.R`, added `R/prone_normalization_steps.R`, `R/zzz_helpers.R`, added `tests/testthat/test-prone_normalization_steps.R`, and added `vignettes/prone_with_probatch.Rmd`.
- Classification/disposition: mixed. Retain only the core QFeatures link guard; exclude PRONE adapters, registrations, tests, vignette, and companion `plot_intragroup_variation` hunk.
- Core hunk/function: `.pb_add_assay_with_link` creates a one-to-one link only when both feature axes are non-null, unique, and set-equal. `pb_transform` and `pb_eval` consume the behavior but are not modified.
- Focused tests/dependencies: companion PRONE chaining/registration tests do not migrate; add a focused core duplicate-row/link test. PRONE moves from `Imports` to `Suggests`, but core must remove it entirely.
- Split targets: `R/ProBatchFeatures.R` for the absent duplicate-row guard; `R/registry.R`, `test-registry.R`, and `test-symbol-ownership.R` for the provider boundary.

### 005 — PRONE condition metadata

<!-- donna-entry kind=commit id=post-sync-003-65f70a46-prone-condition-metadata slot=005 mode=reference-only workflow=- sha=65f70a46c4cf44e2717744aeafbf8acbe83b0378 -->

- Parent/date/message: parent `db128668458a58ad31f66be5b6e39e2fedadbbe1`; `2026-03-14T23:49:32+01:00`; subject `adding condition metadata handling and updating vignette examples to PRONE implementation`; body empty.
- Non-generated paths: `R/impute_PRONE.R`, `R/prone_normalization_steps.R`, `tests/testthat/test-prone_normalization_steps.R`, `vignettes/prone_with_probatch.Rmd`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit contains no core-owned editable artifact or behavior.
- Execution: preserve reserved prefix `005_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks: unifies condition/batch keys, supports combination backends and aliases, validates EigenMS condition requirements, and updates PRONE examples. It affects `imputePRONE*` and PRONE registered normalization providers, not a core export.
- Focused tests/dependencies: PRONE equivalence, chaining, metadata, and backend-compatibility tests; optional `PRONE`, `SummarizedExperiment`, and `S4Vectors`.
- Split targets: provider-aware `R/registry.R`, `test-registry.R`, and `test-symbol-ownership.R`; no provider implementation belongs in core.

### 006 — Grouped missingness heatmaps

<!-- donna-entry kind=commit id=post-sync-004-4e4e9811-grouped-na-heatmaps slot=006 mode=workflow workflow=@/workflows/006_4e4e9811503c_grouped_na_heatmaps.donna.md sha=4e4e9811503cc7da1e35a8113a2fb8383fc5007b -->

- Parent/date/message: parent `65f70a46c4cf44e2717744aeafbf8acbe83b0378`; `2026-03-24T15:26:50+01:00`; subject `Enhance missing-value heatmap functionality with grouped visualization support and improved sample annotation handling`; body empty.
- Non-generated paths: `R/plot_missing.R`, `tests/testthat/test-plot_missing.R`.
- Classification/disposition: core, included.
- Hunks/functions: extracts heatmap, annotation-alignment, and arrangement helpers; adds `plot_grouped_NA_heatmap` generic/default/`ProBatchFeatures` methods; groups by metadata and shows observed fractions and group sizes. Also affects `plot_NA_heatmap.default` and `.ProBatchFeatures`.
- Focused tests/dependencies: one-/multi-column aggregation, default validation, and multi-assay arrangement. Existing `pheatmap`, `gridExtra`, `grid`, `grDevices`, and `SummarizedExperiment`.
- Split targets: `R/plot_missing.R` retains the grouped API; restore the exact grouped tests missing from `test-plot_missing.R`.

### 007 — Split-violin qualification and variancePartition floor

<!-- donna-entry kind=commit id=post-sync-005-12d4597e-violin-qualification slot=007 mode=workflow workflow=@/workflows/007_12d4597e9e86_violin_qualification.donna.md sha=12d4597e9e86ea2daddc77c4d236e692bbfa8fa9 -->

- Parent/date/message: parent `4e4e9811503cc7da1e35a8113a2fb8383fc5007b`; `2026-03-24T15:35:57+01:00`; subject `Update DESCRIPTION and NAMESPACE for variancePartition dependency and add grouped NA heatmap plotting methods`; body empty.
- Non-generated paths: `DESCRIPTION`, `R/plot_split_violinplot.R`.
- Classification/disposition: mixed. Include qualified core calls; exclude companion dependency metadata.
- Core hunks/functions: qualify `ggplot2`, `scales`, and `grid` calls in `GeomSplitViolin`, `geom_split_violin`, and `plot_split_violin_with_boxplot`.
- Companion hunk: raises `variancePartition` to `>= 1.40.1`; that dependency belongs with companion variance-partition exports.
- Focused tests/dependencies: use `test-plot_split_violinplot.R`; no historical test changed. Keep explicit core qualifications and omit `variancePartition`.
- Split targets: `R/plot_split_violinplot.R`, its test, and `DESCRIPTION`; several calls remain unqualified in the pinned split.

### 008 — Grouped heatmap binarization tests

<!-- donna-entry kind=commit id=post-sync-006-5747090a-test-grouped-na-binarization slot=008 mode=workflow workflow=@/workflows/008_5747090a1de1_test_grouped_na_binarization.donna.md sha=5747090a1de1cdde780c6db7a84f988aedb9e8af -->

- Parent/date/message: parent `12d4597e9e86ea2daddc77c4d236e692bbfa8fa9`; `2026-03-24T18:08:32+01:00`; subject `Add tests for .pb_group_missing_matrix binarization and validate force_binarization in plot_grouped_NA_heatmap`; body empty.
- Non-generated path: `tests/testthat/test-plot_missing.R`.
- Classification/disposition: core test-only, included in source order even though commit 9 supplies the implementation.
- Hunks/functions: test `TRUE` as threshold `0.5`, numeric thresholding before `drop_complete`, and rejection outside `[0,1]` for `.pb_group_missing_matrix`/`plot_grouped_NA_heatmap.default`.
- Dependencies: none.
- Split targets: `R/plot_missing.R` already has `force_binarization`; add the absent exact tests after accounting for the next source commit.

### 009 — Grouped density and heatmap binarization

<!-- donna-entry kind=commit id=post-sync-007-7feb6d5c-grouped-na-density slot=009 mode=workflow workflow=@/workflows/009_7feb6d5cea9e_grouped_na_density.donna.md sha=7feb6d5cea9ec92d378453ff65d672f231fa44b8 -->

- Parent/date/message: parent `5747090a1de1cdde780c6db7a84f988aedb9e8af`; `2026-03-24T18:08:49+01:00`; subject `Enhance density plotting functions to support grouped densities and add force_binarization parameter for heatmap visualizations`; body empty.
- Non-generated paths: `R/plot_missing.R`, `tests/testthat/test-plot_missing.R`.
- Classification/disposition: core, included.
- Hunks/functions: implements grouped heatmap binarization and scale selection; extends `plot_NA_density.default`/`.ProBatchFeatures` with metadata grouping, colors, linetypes, parameter forwarding, and faceted multi-assay output.
- Focused tests/dependencies: grouped density data/legends and multi-assay facets plus commit 8 binarization cases; existing `ggplot2`, `grDevices`, `pheatmap`, and `SummarizedExperiment`.
- Split targets: `R/plot_missing.R`; grouped heatmap behavior remains, but public grouped-density parameters and end-to-end tests are absent.

### 010 — Grouped density color schemes

<!-- donna-entry kind=commit id=post-sync-008-0dbce812-density-color-scheme slot=010 mode=workflow workflow=@/workflows/010_0dbce8129444_density_color_scheme.donna.md sha=0dbce8129444c9079a560fb9bcd60b328b1f054c -->

- Parent/date/message: parent `7feb6d5cea9ec92d378453ff65d672f231fa44b8`; `2026-03-24T18:33:29+01:00`; subject `Enhance plot_NA_density functions with color_scheme parameter for grouped densities and update documentation accordingly`; body empty.
- Non-generated paths: `R/plot_missing.R`, `tests/testthat/test-plot_missing.R`.
- Classification/disposition: core, included.
- Hunks/functions: adds brewer/named/list color schemes, raw group keys, common theme, and complete Roxygen parameters for `plot_NA_density.default`, `.ProBatchFeatures`, and helpers.
- Focused tests/dependencies: raw key expectations and named-color verification through `ggplot_build`; existing `ggplot2` and core color helpers.
- Split targets: internal helpers remain in `R/plot_missing.R`, but the public grouped-density signature and test are absent.

### 011 — Sample subsetting and Bayesian PCA

<!-- donna-entry kind=commit id=post-sync-009-6d3c15e4-subset-and-bpca slot=011 mode=workflow workflow=@/workflows/011_6d3c15e41905_subset_and_bpca.donna.md sha=6d3c15e419058d264b61ac43c27e95e0d635ca01 -->

- Parent/date/message: parent `0dbce8129444c9079a560fb9bcd60b328b1f054c`; `2026-03-26T16:13:56+01:00`; subject `Add sample subsetting function and enhance PCA methods with Bayesian support`; body empty.
- Non-generated paths: `DESCRIPTION`, `R/ProBatchFeatures.R`, `R/plot_helpers.R`, `R/proteome_wide_diagnostics.R`, `tests/testthat/test-ProBatchFeatures.R`, `test-plot_helpers.R`, `test-proteome_wide_diagnostics.R`.
- Classification/disposition: core, included with an explicit dependency/API decision for optional Bayesian PCA.
- Hunks/functions: fixes the four-index `ProBatchFeatures` `[` method, adds `pb_subset_samples`, adds per-assay atomic argument handling, and abstracts SVD/BPCA computation for `plot_PCA.default`.
- Focused tests/dependencies: subsetting and assay preservation, atomic title splitting, BPCA with missing data. Adds optional `pcaMethods` to `Suggests`.
- Split targets: `R/ProBatchFeatures.R` already contains `pb_subset_samples` and the class fix; `R/plot_helpers.R` lacks atomic-vector mode; `plot_PCA.default` is SVD-only and `pcaMethods` is absent. Preserve that policy choice for review.

### 012 — omicsGMF all-NA guard and colData compatibility

<!-- donna-entry kind=commit id=post-sync-010-7db05fae-omicsgmf-all-na-guard slot=012 mode=workflow workflow=@/workflows/012_7db05faed18f_omicsgmf_all_na_guard.donna.md sha=7db05faed18f75572d92553ae668288c326bda2a -->

- Parent/date/message: parent `6d3c15e419058d264b61ac43c27e95e0d635ca01`; `2026-03-27T11:45:54+01:00`; subject `Enhance omicsGMF imputation to handle all-NA rows/columns by removing them before fitting and re-inserting as NA in the output`; body empty.
- Non-generated paths: `R/ProBatchFeatures.R`, `R/impute_omicsGMF.R`, `tests/testthat/test-ProBatchFeatures.R`, `test-impute_PRONE.R`, `test-impute_omicsGMF.R`, `test-prone_normalization_steps.R`.
- Classification/disposition: mixed. Include core colData compatibility; exclude provider all-NA logic and provider tests.
- Core hunk/function: `.pb_harmonize_colData` accepts integer-like factor levels equivalent to numeric/integer values while rejecting genuinely different values.
- Companion hunks: omicsGMF removes/restores all-NA axes; PRONE tests ensure no newly all-NA axes. These are reached through companion exports only.
- Focused tests/dependencies: add core factor/integer equivalence and mismatch tests. Companion uses omicsGMF, SingleCellExperiment, SummarizedExperiment, S4Vectors, and PRONE.
- Split targets: `R/ProBatchFeatures.R` lacks the core branch/test; `R/registry.R` and ownership tests enforce exclusion of provider implementations.

### 013 — Centralized all-NA handling

<!-- donna-entry kind=commit id=post-sync-011-fdafc6c8-centralize-all-na-handling slot=013 mode=workflow workflow=@/workflows/013_fdafc6c8fe13_centralize_all_na_handling.donna.md sha=fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5 -->

- Parent/date/message: parent `7db05faed18f75572d92553ae668288c326bda2a`; `2026-03-27T13:35:12+01:00`; subject `Refactor missing value handling: Introduce .pb_strip_allna and .pb_restore_allna functions to streamline removal and restoration of all-NA rows/columns across multiple imputation methods.`; body empty.
- Non-generated paths: `R/ProBatchFeatures.R`, `R/correct_batch_RUVIIIC.R`, `R/correct_batch_effects.R`, `R/fit_non_linear.R`, `R/impute_missForest.R`, `R/impute_omicsGMF.R`, `R/normalize.R`, `R/utility_funcs.R`.
- Classification/disposition: mixed. Review each core hunk; exclude RUVIIIC, missForest, and omicsGMF consumers.
- Core hunks/functions: `.pb_strip_allna`/`.pb_restore_allna`; no-handling boolean semantics for registry median normalization and correction; wide validation in `center_feature_batch`; deprecated forwarding for `correct_batch_effects_df/dm`; nullable group handling in `normalize_sample_medians_df`; comment-only `fit_nonlinear` clarification.
- Focused tests/dependencies: none historically. Compare `test-correct_batch_effect.R`, `test-batch_effect_steps.R`, `test-handle_missing_values.R`, `test-matrix_adapter.R`, `test-normalize.R`, `test-fit_non_linear.R`, and `test-registry.R`. Provider dependencies stay out of core.
- Split targets: wide validation is present; canonical missing policy supersedes booleans; all-NA helpers are absent; group-aware median normalization was removed; deprecated wrapper behavior needs explicit API review.

### 014 — Demote all-NA helper Roxygen

<!-- donna-entry kind=commit id=post-sync-012-77a3365f-demote-all-na-helper-roxygen slot=014 mode=workflow workflow=@/workflows/014_77a3365f0e13_demote_all_na_helper_roxygen.donna.md sha=77a3365f0e13e9abca40a14e1898a351a1cdbc7f -->

- Parent/date/message: parent `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5`; `2026-03-27T13:38:21+01:00`; subject `Refactor documentation for .pb_strip_allna and .pb_restore_allna functions to improve clarity and consistency`; body empty.
- Non-generated path: `R/utility_funcs.R`.
- Classification/disposition: metadata/documentation-only. Apply only if the preceding helper design is adopted.
- Hunk: converts the two internal helper Roxygen blocks to ordinary comments; no runtime or exported API change.
- Tests/dependencies: none.
- Split targets: helpers are absent, so no current target exists.

### 015 — Group masking and NA-intensity plot

<!-- donna-entry kind=commit id=post-sync-013-4b117a03-mask-groups-plot-na-intensity slot=015 mode=workflow workflow=@/workflows/015_4b117a03f0c5_mask_groups_plot_na_intensity.donna.md sha=4b117a03f0c5295f2253db26e13b945d7e3027b5 -->

- Parent/date/message: parent `77a3365f0e13e9abca40a14e1898a351a1cdbc7f`; `2026-03-28T12:00:50+01:00`; subject `Add mask_failing parameter to pb_groupfilterNA and implement plot_NA_intensity function for visualizing missingness versus mean intensity`; body empty.
- Non-generated paths: `R/pb_missing_filters.R`, added `R/plot_NA_intensity.R`, `tests/testthat/test-pb_missing_helpers.R`, added `tests/testthat/test-plot_NA_intensity.R`.
- Classification/disposition: core, included.
- Hunks/functions: adds/logs `mask_failing`; adds `plot_NA_intensity` generic/default/`ProBatchFeatures` methods, grouped statistics, spline trend, Spearman labels, and facets.
- Focused tests/dependencies: grouped masking and NA-intensity tests; existing SummarizedExperiment, ggplot2, scales, and splines.
- Split targets: both implementations exist in `R/pb_missing_filters.R` and `R/plot_NA_intensity.R`; restore two missing focused `mask_failing` cases.

### 016 — Clarify removeBatchEffect missing policy

<!-- donna-entry kind=commit id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy slot=016 mode=workflow workflow=@/workflows/016_4285c42f3167_clarify_removebatcheffect_missing_policy.donna.md sha=4285c42f31670d2f750dc8eb8c7ff1d0134a342d -->

- Parent/date/message: parent `4b117a03f0c5295f2253db26e13b945d7e3027b5`; `2026-03-31T14:20:27+02:00`; subject `Refactor missing value handling: Update warning messages and improve NA handling logic in correct_with_removeBatchEffect_dm and handle_missing_values functions.`; body empty.
- Non-generated paths: `R/correct_batch_effects.R`, `R/handle_missing_values.R`.
- Classification/disposition: core, include semantics rather than legacy boolean syntax.
- Hunks/functions: improves warnings, skips preprocessing when false, emits removal warnings only when removals occur, and removes an inaccurate limma claim. Affects `correct_with_removeBatchEffect_dm`, `.removeBatchEffect_matrix_step`, `.run_matrix_method`, and `handle_missing_values`.
- Focused tests/dependencies: no historical test; compare missing/correction/matrix-adapter tests. No dependency change.
- Split targets: `R/handle_missing_values.R`, `R/matrix_adapter.R`, and `R/correct_batch_effects.R`; canonical `error/keep/drop_features/fill` plus `fill_value` supersedes the boolean interface.

### 017 — Generated missing-value docs

<!-- donna-entry kind=commit id=post-sync-015-4540aca9-generated-missing-docs slot=017 mode=workflow workflow=@/workflows/017_4540aca9182c_generated_missing_docs.donna.md sha=4540aca9182c6708fe9bda0b8fc33d2cf8c13e57 -->

- Parent/date/message: parent `4285c42f31670d2f750dc8eb8c7ff1d0134a342d`; `2026-03-31T14:20:33+02:00`; subject `Refactor documentation for handle_missing_values function: Improve clarity by removing redundant mention of limma's removeBatchEffect.`; body empty.
- Non-generated paths: none.
- Classification/disposition: generated-only, explicit no-change child. Do not inspect or migrate generated output.
- Core functions/tests/dependencies: none beyond the editable-source change already accounted for in commit 16.
- Split target: none; maintainers regenerate documentation later.

### 018 — omicsGMF batch-mean semantics

<!-- donna-entry kind=commit id=post-sync-016-b8a262b4-subtract-omicsgmf-batch-mean slot=018 mode=reference-only workflow=- sha=b8a262b4256966d60e1e8452ebde7a1bf471b4af -->

- Parent/date/message: parent `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57`; `2026-05-22T11:54:56+02:00`; subject `Refactor batch correction logic in correct_with_omicsGMF: Update to subtract only the batch-attributable mean from observed data, aligning with limma::removeBatchEffect semantics. Adjust tests to reflect new model semantics for non-batch design terms.`; body empty.
- Non-generated paths: `R/correct_batch_omicsGMF.R`, `tests/testthat/test-correct_batch_omicsGMF.R`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit contains no core-owned editable artifact or behavior.
- Execution: preserve reserved prefix `018_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks: `correct_with_omicsGMF` subtracts only batch design contribution, validates Beta orientation/dimensions, and preserves latent-only fallback.
- Focused tests/dependencies: observed-minus-batch-mean semantics; omicsGMF/sgdGMF/SingleCellExperiment stack.
- Split targets: only `R/registry.R`, `R/matrix_adapter.R`, and ownership tests; implementation belongs in Bench.

### 019 — CV, dates, LOESS, and mComBat safeguards

<!-- donna-entry kind=commit id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat slot=019 mode=workflow workflow=@/workflows/019_d95dd736cb27_harden_cv_dates_loess_mcombat.donna.md sha=d95dd736cb27d68fb2e21b20ab97d0f69826a663 -->

- Parent/date/message: parent `b8a262b4256966d60e1e8452ebde7a1bf471b4af`; `2026-05-22T11:55:44+02:00`; subject `Enhance CV calculation and m-ComBat error handling: Implement safeguards against division by zero in compute_cv and enforce minimum sample size per batch in correct_with_mComBat. Add warnings for features with zero variance and suppress extrapolation in LOESS regression functions.`; body empty.
- Non-generated paths: `R/CV_calculation.R`, `R/M_ComBat.R`, `R/date_conversion.R`, `R/fit_non_linear.R`.
- Classification/disposition: mixed. Include all four core safeguards; exclude mComBat.
- Core hunks/functions: `compute_cv` returns `NA` for non-finite/near-zero means; `dates_to_posix` restores locale; `date_to_sample_order` uses first-tie order; LOESS helpers suppress extrapolation and return NA fallbacks.
- Companion hunk: mComBat validates batch sizes and within-batch variance and leaves invalid features unchanged with warnings.
- Focused tests/dependencies: add cases to `test-CV_calculation.R`, `test-date_conversion.R`, and `test-fit_non_linear.R`; no core dependency change.
- Split targets: all four core safeguards are absent in the corresponding pinned sources. mComBat is absent intentionally.

### 020 — Sample alignment and consistency

<!-- donna-entry kind=commit id=post-sync-018-6989542e-harden-sample-alignment-consistency slot=020 mode=workflow workflow=@/workflows/020_6989542e4ece_harden_sample_alignment_consistency.donna.md sha=6989542e4ece84f9fe210f3bf06b8254684e425e -->

- Parent/date/message: parent `d95dd736cb27d68fb2e21b20ab97d0f69826a663`; `2026-05-22T14:37:21+02:00`; subject `Enhance sample consistency checks and warnings: Improve feedback on mismatches between sample annotation and data matrix, detailing samples unique to each. Update corrected matrix handling to warn on sample set differences and reorder surviving columns accordingly.`; body empty.
- Non-generated paths: `R/correct_batch_effects.R`, `R/correlation-based_diagnostics.R`, `R/fit_non_linear.R`, `R/impute_PRONE.R`, `R/normalize.R`, `R/utility_funcs.R`.
- Classification/disposition: mixed. Include independent core hunks; exclude PRONE and deliberately removed group-aware normalization.
- Core hunks/functions: detailed correction sample warnings/order restoration; pairwise `plot_sample_corr_heatmap`; caller-copy protection in `fit_nonlinear`; skipped-group reporting in group-aware normalization; detailed `check_sample_consistency`.
- Focused tests/dependencies: correction, correlation, fit, and utility tests; normalization only if its retained API needs it. No dependency change.
- Split targets: correction reordering is already present; pairwise correlation, caller-copy protection, and detailed consistency warning are absent. PRONE has no core target.

### 021 — Provider control-flow and source-doc cleanup

<!-- donna-entry kind=commit id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs slot=021 mode=workflow workflow=@/workflows/021_428ad74b73fb_cleanup_provider_control_flow_docs.donna.md sha=428ad74b73fbc067de09d0a71c93395ae85eb51f -->

- Parent/date/subject: parent `6989542e4ece84f9fe210f3bf06b8254684e425e`; `2026-05-22T20:03:22+02:00`; complete subject `Refactor various functions and documentation: - Replace for loops with seq_len for better performance in correct_with_mComBat and correct_with_PLSDA_batch. - Update warning messages and improve NA handling in correct_with_NormAE. - Enhance documentation with return values for several functions. - Modify examples to use \donttest for non-executable code snippets in plot functions.`; body empty.
- Non-generated paths: `R/M_ComBat.R`, `R/ProBatchFeatures.R`, `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_NormAE.R`, `R/correct_batch_RUVIIIC.R`, `R/correct_batch_effects_old.R`, `R/fit_non_linear.R`, `R/impute_PRONE.R`, `R/plot_NA_intensity.R`, `R/proBatch.R`, `R/proteome_wide_diagnostics.R`.
- Classification/disposition: mixed. Include core Roxygen/comment/control-flow hunks; exclude every provider and TSNE/UMAP hunk.
- Core hunks/functions: return docs for the class/package/legacy wrappers, environment-based LOESS warning capture, and a `\donttest` example for `plot_NA_intensity`.
- Companion hunks: mComBat, PLSDA, NormAE, RUVIIIC, PRONE, TSNE, and UMAP cleanup.
- Focused tests/dependencies: no historical tests; compare class, correction, fit, and NA-intensity tests. No metadata change.
- Split targets: NA-intensity already uses `\donttest`; LOESS capture and some return docs are absent; correction docs moved to consolidated `R/correct_batch_effects.R`.

### 022 — Version and variancePartition metadata

<!-- donna-entry kind=commit id=post-sync-020-72d11d1f-version-variancepartition slot=022 mode=workflow workflow=@/workflows/022_72d11d1f7f92_version_variancepartition.donna.md sha=72d11d1f7f92a10cb9a5244e68f9c19e830197b9 -->

- Parent/date/message: parent `428ad74b73fbc067de09d0a71c93395ae85eb51f`; `2026-05-22T22:38:23+02:00`; subject `Bump version to 2.1.0 and update variancePartition dependency to >= 1.36.2`; body empty.
- Non-generated path: `DESCRIPTION`.
- Classification/disposition: metadata. Do not mechanically port either hunk.
- Hunks: version `1.99.5` to `2.1.0`; actually relaxes `variancePartition` from `>= 1.40.1` to `>= 1.36.2`.
- Core functions/tests: none; release metadata review only.
- Split targets: current destination already has version `2.1.0`; `variancePartition` is companion-owned and absent from the pinned split core.

### 023 — Docker ignore metadata

<!-- donna-entry kind=commit id=post-sync-021-8088de17-ignore-docker-readme slot=023 mode=workflow workflow=@/workflows/023_8088de1701e0_ignore_docker_readme.donna.md sha=8088de1701e0908cca25b978105f8d6c7bfccc20 -->

- Parent/date/message: parent `72d11d1f7f92a10cb9a5244e68f9c19e830197b9`; `2026-05-23T20:36:23+02:00`; subject `Update .gitignore to include Docker-related files and directories`; body empty.
- Non-generated path: `.gitignore`.
- Classification/disposition: metadata, explicit review/no-change candidate.
- Hunk: reorders existing Docker ignores and adds `Docker_README.md`; no API, test, or dependency effect.
- Split target: pinned `.gitignore` lacks `Docker_README.md` but adds `Rplots.pdf`; destination policy controls the result.

### 024 — Companion missing-policy forwarding

<!-- donna-entry kind=commit id=post-sync-022-4a06d999-forward-companion-missing-policy slot=024 mode=reference-only workflow=- sha=4a06d99949114b2804b8d34492a288872fb611ed -->

- Parent/date/subject: parent `8088de1701e0908cca25b978105f8d6c7bfccc20`; `2026-05-23T20:37:10+02:00`; subject `Enhance handling of missing values across multiple functions`.
- Complete body:
  - `Added fill_the_missing parameter to correct_with_omicsGMF, imputePRONE_dm, missForestImpute, and estimate_omicsGMF_rank functions to allow for more flexible handling of missing data.`
  - `Updated internal methods to forward fill_the_missing parameter appropriately, ensuring consistent behavior across different imputation and correction methods.`
  - `Implemented regression tests to verify that fill_the_missing is correctly passed through and utilized in matrix steps for omicsGMF, PRONE, and missForest.`
  - `Adjusted documentation for relevant functions to reflect the new parameter and its usage.`
- Non-generated paths: `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_omicsGMF.R`, `R/impute_PRONE.R`, `R/impute_missForest.R`, `R/impute_omicsGMF.R`, `R/prone_normalization_steps.R`, and matching `test-correct_batch_BERT_PLSDA.R`, `test-correct_batch_omicsGMF.R`, `test-impute_PRONE.R`, `test-impute_missForest.R`, `test-impute_omicsGMF.R`, and `test-prone_normalization_steps.R`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit changes only provider adapters and their tests, not core missing-value functions.
- Execution: preserve reserved prefix `024_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks/functions: forwards `fill_the_missing` through BERT, PLSDA, omicsGMF correction/imputation/rank, PRONE imputation/normalization, and missForest paths.
- Focused tests/dependencies: forwarding and positional-binding cases across those providers; BERT, BiocParallel, PLSDAbatch, mixOmics, PRONE, missForest, omicsGMF, sgdGMF, SingleCellExperiment.
- Split targets: translate only the boundary concept to canonical `missing`/`fill_value` in `R/matrix_adapter.R`, `R/handle_missing_values.R`, `R/registry.R`, and `R/step_result.R`; enforce negative ownership tests.

### 025 — omicsGMF design names

<!-- donna-entry kind=commit id=post-sync-023-6601232c-restore-omicsgmf-design-names slot=025 mode=reference-only workflow=- sha=6601232c69b44c507f2e9f63d256836e655f7973 -->

- Parent/date/message: parent `4a06d99949114b2804b8d34492a288872fb611ed`; `2026-05-23T21:56:56+02:00`; subject `Fix column name restoration in design matrix for omicsGMF fitting`; body empty.
- Non-generated path: `R/impute_omicsGMF.R`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit contains no core-owned editable artifact or behavior.
- Execution: preserve reserved prefix `025_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunk: restores design-matrix column names after `sgdGMF` strips them so named batch terms can be matched.
- Tests/dependencies: no test added, an explicit coverage gap; omicsGMF, sgdGMF, SingleCellExperiment.
- Split targets: no provider implementation; boundary-only comparison with `R/matrix_adapter.R`, `R/step_result.R`, and ownership tests.

### 026 — omicsGMF latent artifacts

<!-- donna-entry kind=commit id=post-sync-024-96d38eb7-preserve-omicsgmf-latents slot=026 mode=reference-only workflow=- sha=96d38eb7449e4d38c0d6a3fffd66e17f145f669f -->

- Parent/date/message: parent `6601232c69b44c507f2e9f63d256836e655f7973`; `2026-05-24T16:55:41+02:00`; subject `Attach omicsGMF latent representation as attributes to corrected and imputed matrices for downstream analysis`; body empty.
- Non-generated paths: `R/correct_batch_omicsGMF.R`, `R/impute_omicsGMF.R`, `test-correct_batch_omicsGMF.R`, `test-impute_omicsGMF.R`.
- Classification/disposition: companion-only, reference-only skip. The source implementation does not enter core; its structured-artifact implication remains covered by the residual split review.
- Execution: preserve reserved prefix `026_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks/functions: attach scores, loadings, design `X`, design `Beta`, and reduction name to correction/imputation matrices.
- Focused tests/dependencies: correction/imputation attribute tests; omicsGMF, sgdGMF, SingleCellExperiment, SummarizedExperiment, S4Vectors.
- Split targets: `R/step_result.R::pb_step_result`, `R/ProBatchFeatures.R::pb_transform`, and `test-step_result.R`; structured artifacts supersede matrix attributes.

### 027 — Explicit final log assay

<!-- donna-entry kind=commit id=post-sync-025-b60e2b16-store-explicit-final-log-assay slot=027 mode=workflow workflow=@/workflows/027_b60e2b169afe_store_explicit_final_log_assay.donna.md sha=b60e2b169afe45934c18eee1c83469e7e5fed33f -->

- Parent/date/message: parent `96d38eb7449e4d38c0d6a3fffd66e17f145f669f`; `2026-05-24T22:40:21+02:00`; subject `Enhance log transformation functions and update pb_transform to conditionally store transformed values with final_name`; body empty.
- Non-generated paths: `R/ProBatchFeatures.R`, `tests/testthat/test-ProBatchFeatures.R`.
- Classification/disposition: core, include or record already-equivalent behavior.
- Hunks/functions: log/log2 registry steps accept `log_base`, `offset`, and legacy aliases; `pb_transform` materializes an otherwise ephemeral final log step when `final_name` is supplied without renaming the source assay.
- Focused tests/dependencies: transformed-not-raw final assay; base R and existing QFeatures/SummarizedExperiment/S4Vectors infrastructure.
- Split targets: `R/registry.R`, `R/ProBatchFeatures.R`, and `test-ProBatchFeatures.R`; behavior is already absorbed and strengthened with name validation and lineage synchronization.

### 028 — In-place filter link effects

<!-- donna-entry kind=commit id=post-sync-026-5edcc4b8-report-inplace-link-effects slot=028 mode=workflow workflow=@/workflows/028_5edcc4b89f5c_report_inplace_link_effects.donna.md sha=5edcc4b89f5cbfba07e1ca1057ab475a56501202 -->

- Parent/date/message: parent `b60e2b169afe45934c18eee1c83469e7e5fed33f`; `2026-05-24T23:09:39+02:00`; subject `Add warning for in-place filtering in pb_filterNA and pb_groupfilterNA to inform about QFeatures assay link removal`; body empty.
- Non-generated path: `R/pb_missing_filters.R`.
- Classification/disposition: core, but split behavior has evolved.
- Hunks/functions: `pb_filterNA`/`pb_groupfilterNA` detect links and warn that in-place filtering removes them; Roxygen explains the effect.
- Focused tests/dependencies: no historical test; QFeatures and SummarizedExperiment.
- Split targets: `R/pb_missing_filters.R` and `test-pb_missing_helpers.R`; split `pb_filterNA` preserves/prunes valid links, while `pb_groupfilterNA` retains a removal warning. Do not replay obsolete warnings.

### 029 — PRONE default input assay

<!-- donna-entry kind=commit id=post-sync-027-0b3c3b55-default-prone-input-to-log2 slot=029 mode=reference-only workflow=- sha=0b3c3b55403f2bb7342a37236dd30f3d9b3544e9 -->

- Parent/date/message: parent `5edcc4b89f5cbfba07e1ca1057ab475a56501202`; `2026-05-25T00:54:40+02:00`; subject `Update assay input parameter from "raw" to "log2" in normalization functions and tests`; body empty.
- Non-generated paths: `R/prone_normalization_steps.R`, `tests/testthat/test-prone_normalization_steps.R`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit contains no core-owned editable artifact or behavior.
- Execution: preserve reserved prefix `029_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks: PRONE registered normalization defaults `assay_in` to `log2`; raw-behavior tests become explicit.
- Dependencies/tests: PRONE, SummarizedExperiment, S4Vectors and provider tests.
- Split targets: provider-aware registry and ownership tests only; no PRONE implementation in core.

### 030 — Combined omicsGMF imputation/correction

<!-- donna-entry kind=commit id=post-sync-028-20e76c9b-combine-omicsgmf-impute-correct slot=030 mode=reference-only workflow=- sha=20e76c9b9b28c0ec98faa63c3f9382c1347301b9 -->

- Parent/date/message: parent `0b3c3b55403f2bb7342a37236dd30f3d9b3544e9`; `2026-05-25T03:00:22+02:00`; subject `Add impute_and_correct_with_omicsGMF function for combined imputation and batch correction`; body empty.
- Non-generated paths: `R/correct_batch_omicsGMF.R`, `R/impute_omicsGMF.R`, `R/zzz_helpers.R`, `tests/testthat/test-correct_batch_omicsGMF.R`.
- Classification/disposition: companion-only, reference-only skip. No executable child workflow is warranted because the commit contains no core-owned editable artifact or behavior.
- Execution: preserve reserved prefix `030_`, retain this record in source order, never run a child for it, and record only a deterministic `reference-only` progress outcome when reached.
- Hunks/functions: adds exported `impute_and_correct_with_omicsGMF`, reuses one fit, subtracts batch contribution from imputed values, adds fallback control, fixes forwarding, and registers aliases.
- Focused tests/dependencies: delegation, corrected imputed values, fallback semantics; omicsGMF, sgdGMF, SingleCellExperiment, SummarizedExperiment, S4Vectors.
- Split targets: `R/registry.R`, `R/matrix_adapter.R`, `R/step_result.R`, and their tests. Export and aliases belong exclusively in Bench.

### 031 — Unused model levels

<!-- donna-entry kind=commit id=post-sync-029-e2bb1854-drop-unused-model-levels slot=031 mode=workflow workflow=@/workflows/031_e2bb18547c73_drop_unused_model_levels.donna.md sha=e2bb18547c73f1c471fc1afcb3facbd8bea5fa92 -->

- Parent/date/message: parent `20e76c9b9b28c0ec98faa63c3f9382c1347301b9`; `2026-05-25T17:13:44+02:00`; subject `Enhance BERT batch correction functions to handle unused factor levels and expand categorical covariates into dummy variables for improved modeling accuracy`; body empty.
- Non-generated paths: `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_effects.R`, `tests/testthat/test-bert-covariates.R`, `tests/testthat/test-correct_batch_effect.R`.
- Classification/disposition: mixed. Include unused-level handling for core correction; exclude all BERT implementation and tests.
- Core hunks/functions: `run_ComBat_core` and `.removeBatchEffect_matrix_step` drop unused batch/covariate levels before model construction, transitively affecting `correct_with_ComBat*` and `correct_with_removeBatchEffect*`.
- Companion hunks: BERT drops unused batch levels and expands multi-level categorical covariates to `k-1` dummy columns while retaining one-column binary coding.
- Focused tests/dependencies: core unused-factor regression in `test-correct_batch_effect.R`; core `sva`, `limma`, and `stats`; companion BERT/BiocParallel remain outside core.
- Split targets: the pinned core correction implementations and test lack `droplevels` handling, so this core hunk remains actionable. The BERT target is only the registry/ownership boundary.

### 032 — Residual stopped-split review

<!-- donna-entry kind=residual id=residual-split-review slot=032 mode=workflow workflow=@/workflows/032_residual_split_review.donna.md sha=- -->

- Classification/disposition: residual evidence and decision workflow, always included after the foundation, 22 executable commit children, and seven recorded reference-only skips account for all 29 source commits.
- Scope: compare the completed destination with the complete stopped split and provenance note, including the mode-only preflight and every non-generated baseline-to-split artifact family.
- Exact residual inventory: use pinned objects for `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92..29a7478dc7deea846a2c1ff1abd25a881e6f87db`, excluding `man/**` and `NAMESPACE`. The expected immutable inventory is 99 paths: 11 added, 38 deleted, and 50 modified, with name-status SHA-256 `28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38`.
- Required decisions: the split-only identifier, matrix-adapter, registry, and structured-result layers; associated focused tests; `pb_apply_matrix_method`, `pb_step_result`, and `pb_unregister_steps`; all 26 companion exports; provider removals; duplicate-definition consolidation; and retained/deleted coverage.
- Outputs: `workflows/reports/bec-ecoli-core-residual-review.md` and `workflows/reports/bec-ecoli-core-remaining-change-plan.md`.
- Review control: its report-review request remains pending until the maintainer explicitly resumes after manual analysis and any manual report commit.
