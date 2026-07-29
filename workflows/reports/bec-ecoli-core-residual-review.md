# BEC E. coli Core Residual Split Review

## Scope and evidence

This report audits the stopped split only as historical evidence after
workflows `002` through `031` completed. The stopped provenance record identifies
`29a7478dc7deea846a2c1ff1abd25a881e6f87db` as an unaccepted implementation,
not a merge recommendation. The mode-only preflight
`49cee7cc978fbb149c262a5a783face32dd1d135` changed `DESCRIPTION` from `100755`
to `100644`; the current destination already uses the regular-file mode.

The non-generated baseline-to-split comparison
`ba6ee246eace090e71baa7aba302ca64e76ddb32..29a7478dc7deea846a2c1ff1abd25a881e6f87db`
contains 47 paths (19 added and 28 modified; 9,489 insertions and 677
deletions). The exhaustive residual comparison
`e2bb18547c73f1c471fc1afcb3facbd8bea5fa92..29a7478dc7deea846a2c1ff1abd25a881e6f87db`
contains 99 paths: 11 added, 38 deleted, and 50 modified. Its exact
name-status stream has SHA-256
`28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38`.
Generated manuals and `NAMESPACE` are excluded from every comparison.

Classification totals are:

- `required`: 4 paths;
- `recommended`: 7 paths;
- `decision`: 13 paths;
- `equivalent`: 43 paths;
- `excluded`: 32 paths.

`Equivalent` includes destination behavior that is already stronger than the
stopped split. `Excluded` is reserved for companion-owned implementation,
tests, resources, and vignettes that must remain outside core.

## Immutable 99-path inventory

### Project metadata

<!-- split-path status=M path=.Rbuildignore class=equivalent -->
<!-- split-path status=M path=.gitignore class=equivalent -->
<!-- split-path status=M path=DESCRIPTION class=decision -->
<!-- split-path status=M path=NEWS class=required -->

The current ignore files preserve project workflow/session exclusions and
already reject the split's unsupported `Rplots.pdf` policy. `DESCRIPTION`
already excludes provider dependencies and intentionally keeps version `2.1.0`
and R `>= 4.5.0`; the split's `1.99.5` and R `>= 4.6.0` values are not release
authority. A maintainer decision remains on the release-version signal and on
optional dependency cleanup after the focused `reshape2` and `plyr` removals.
`NEWS` must then document the accepted identifier, missing-policy,
group-masking, plotting, CV, LOESS, lineage, and correction changes. Validate
metadata as DCF, run the focused dependency checks, and perform release checks
only after those decisions.

### R implementation

<!-- split-path status=M path=R/CV_calculation.R class=equivalent -->
<!-- split-path status=D path=R/M_ComBat.R class=excluded -->
<!-- split-path status=M path=R/ProBatchFeatures.R class=required -->
<!-- split-path status=M path=R/auxiliary.R class=equivalent -->
<!-- split-path status=D path=R/batch_correction_helpers.R class=recommended -->
<!-- split-path status=D path=R/calculate_intragroup.R class=excluded -->
<!-- split-path status=D path=R/classification_metrics.R class=excluded -->
<!-- split-path status=M path=R/colors_for_annotation.R class=equivalent -->
<!-- split-path status=D path=R/correct_batch_BERT_PLSDA.R class=excluded -->
<!-- split-path status=D path=R/correct_batch_NormAE.R class=excluded -->
<!-- split-path status=D path=R/correct_batch_RUVIIIC.R class=excluded -->
<!-- split-path status=M path=R/correct_batch_effects.R class=decision -->
<!-- split-path status=D path=R/correct_batch_effects_old.R class=equivalent -->
<!-- split-path status=D path=R/correct_batch_omicsGMF.R class=excluded -->
<!-- split-path status=M path=R/correlation-based_diagnostics.R class=equivalent -->
<!-- split-path status=M path=R/date_conversion.R class=equivalent -->
<!-- split-path status=D path=R/explained_variance_plots.R class=equivalent -->
<!-- split-path status=M path=R/feature_level_diagnostics.R class=equivalent -->
<!-- split-path status=M path=R/fit_non_linear.R class=equivalent -->
<!-- split-path status=M path=R/handle_missing_values.R class=equivalent -->
<!-- split-path status=A path=R/identifiers.R class=decision -->
<!-- split-path status=D path=R/impute_PRONE.R class=excluded -->
<!-- split-path status=D path=R/impute_missForest.R class=excluded -->
<!-- split-path status=D path=R/impute_omicsGMF.R class=excluded -->
<!-- split-path status=M path=R/initial_assessment.R class=equivalent -->
<!-- split-path status=A path=R/matrix_adapter.R class=decision -->
<!-- split-path status=M path=R/metadata_diagnostics.R class=equivalent -->
<!-- split-path status=M path=R/normalize.R class=equivalent -->
<!-- split-path status=M path=R/pb_missing_filters.R class=equivalent -->
<!-- split-path status=M path=R/pbf_input_helpers.R class=recommended -->
<!-- split-path status=M path=R/plot_NA_intensity.R class=equivalent -->
<!-- split-path status=M path=R/plot_helpers.R class=equivalent -->
<!-- split-path status=M path=R/plot_missing.R class=equivalent -->
<!-- split-path status=M path=R/plot_split_violinplot.R class=recommended -->
<!-- split-path status=M path=R/proBatch.R class=decision -->
<!-- split-path status=D path=R/prone_normalization_steps.R class=excluded -->
<!-- split-path status=M path=R/proteome_wide_diagnostics.R class=equivalent -->
<!-- split-path status=A path=R/registry.R class=decision -->
<!-- split-path status=A path=R/step_result.R class=decision -->
<!-- split-path status=M path=R/transform_raw.R class=equivalent -->
<!-- split-path status=M path=R/utility_funcs.R class=equivalent -->
<!-- split-path status=M path=R/zzz_helpers.R class=decision -->

The ordered migration already adopted or superseded the stopped split's core
CV, conversion, color, correlation, date, diagnostic, missing-value,
normalization, plotting, transformation, PVCA, and utility behavior. The
current group-aware median normalization and qualified plotting calls are
stronger than the comparator's reduced forms. Core PVCA ownership is already
consolidated in `R/proteome_wide_diagnostics.R`, and the obsolete
`R/correct_batch_effects_old.R` file is absent.

Two concrete lineage defects remain in `R/ProBatchFeatures.R`:
`pb_add_level()` silently accepts an existing target without checking the
requested data and parent lineage, and `.pb_add_log_entry()` permits a second,
conflicting non-self parent edge for an existing result. These are independent
of the rejected registry design and require focused conflict checks first.

`R/correct_batch_effects.R` still defines `correct_batch_effects_df`,
`correct_batch_effects_dm`, and `correct_with_removeBatchEffect_dm` twice. The
effective later definitions are known, but deleting the shadowed copies and
adding a repository-wide uniqueness invariant remains a maintainer decision.
`R/batch_correction_helpers.R` is intentionally retained as a focused core
helper; its remaining action is a Depmesh coverage relation, not split-style
inlining. The same governance recommendation applies to
`R/pbf_input_helpers.R`.

The current identifier hardening is accepted. The comparator's additional
positional-rowname prohibition, mandatory feature IDs, and opt-in duplicate
aggregation remain a coupled API decision. The absent matrix adapter,
provider-aware registry, and structured-step result add the public
`pb_apply_matrix_method()`, `pb_unregister_steps()`, and `pb_step_result()`
contracts. Decide them in that dependency order and preserve the accepted
`loessLimmaRBE` registrations if the registry is redesigned. `R/zzz_helpers.R`
changes only with that provider-lifecycle decision.

The split demonstrates that `plyr::arrange()` can be replaced by base ordering
inside `R/plot_split_violinplot.R`, but its de-qualified ggplot calls are
rejected. `R/proBatch.R` remains a decision because dependency cleanup must
keep every still-used import and its Roxygen source synchronized with
`DESCRIPTION`.

The deleted mComBat, BERT, NormAE, RUVIIIC, omicsGMF, PRONE, missForest,
classification, intragroup, and provider-normalization sources are explicitly
companion-owned. Their dependencies and implementation must not re-enter core.

### User-facing and installed resources

<!-- split-path status=M path=README.md class=equivalent -->
<!-- split-path status=M path=inst/CITATION class=recommended -->
<!-- split-path status=D path=inst/overrides/PRONE/Imputation.R class=excluded -->
<!-- split-path status=D path=inst/overrides/normae/__main__.py class=excluded -->
<!-- split-path status=D path=inst/scripts/README.md class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_task_helpers.R class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_workflow.R class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_workflow_params.R class=excluded -->
<!-- split-path status=D path=inst/scripts/pb_tasks.yaml class=excluded -->

`README.md` is byte-identical to the comparator. `inst/CITATION` has an
independent maintenance improvement: replace deprecated `citEntry()` with
`bibentry()` and correct `Evan G. Willams` to `Evan G. Williams`, then parse
the citation file. The vendored PRONE override, NormAE Python overlay, and
benchmark scripts/YAML belong exclusively to the companion workflow and
remain excluded.

### Focused tests

<!-- split-path status=A path=tests/testthat/helper-source-root.R class=recommended -->
<!-- split-path status=M path=tests/testthat/test-CV_calculation.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-M_ComBat.R class=excluded -->
<!-- split-path status=M path=tests/testthat/test-ProBatchFeatures.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-ProBatchFeatures_links.R class=recommended -->
<!-- split-path status=M path=tests/testthat/test-auxiliary.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-batch_effect_steps.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-bert-covariates.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-calculate_intragroup.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-classification_metrics.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-correct_batch_BERT_PLSDA.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-correct_batch_NormAE.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-correct_batch_RUVIIIC.R class=excluded -->
<!-- split-path status=M path=tests/testthat/test-correct_batch_effect.R class=required -->
<!-- split-path status=D path=tests/testthat/test-correct_batch_omicsGMF.R class=excluded -->
<!-- split-path status=M path=tests/testthat/test-correlation_based_diagnostics.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-design_diagnostics.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-explained_variance_plots.R class=recommended -->
<!-- split-path status=M path=tests/testthat/test-feature_level_diagnostics.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-fit_non_linear.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-handle_missing_values.R class=equivalent -->
<!-- split-path status=A path=tests/testthat/test-identifiers.R class=decision -->
<!-- split-path status=D path=tests/testthat/test-impute_PRONE.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-impute_missForest.R class=excluded -->
<!-- split-path status=D path=tests/testthat/test-impute_omicsGMF.R class=excluded -->
<!-- split-path status=M path=tests/testthat/test-initial_assessment.R class=equivalent -->
<!-- split-path status=A path=tests/testthat/test-lineage.R class=required -->
<!-- split-path status=A path=tests/testthat/test-matrix_adapter.R class=decision -->
<!-- split-path status=M path=tests/testthat/test-metadata_diagnostics.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-normalize.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-pb_missing_helpers.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-plot_NA_intensity.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-plot_helpers.R class=equivalent -->
<!-- split-path status=M path=tests/testthat/test-plot_missing.R class=equivalent -->
<!-- split-path status=D path=tests/testthat/test-prone_normalization_steps.R class=excluded -->
<!-- split-path status=M path=tests/testthat/test-proteome_wide_diagnostics.R class=equivalent -->
<!-- split-path status=A path=tests/testthat/test-registry.R class=decision -->
<!-- split-path status=A path=tests/testthat/test-step_result.R class=decision -->
<!-- split-path status=A path=tests/testthat/test-symbol-ownership.R class=decision -->
<!-- split-path status=M path=tests/testthat/test-utility_funcs.R class=equivalent -->

The accepted migration has focused coverage for the equivalent core behavior;
retained link, PVCA, and plot-helper suites must not be deleted wholesale.
`test-ProBatchFeatures_links.R` and `test-explained_variance_plots.R` instead
need explicit Depmesh mappings. `test-correct_batch_effect.R` must replace its
undeclared `reshape2::dcast()` use with an already-declared or package-owned
conversion.

Selectively port the two conflict cases from `test-lineage.R` into the current
container suite: existing result data/lineage mismatch and multiple parent
edges. Other lineage cases are already covered after workflow `027`.
`helper-source-root.R` is useful shared infrastructure for source-ownership
tests.

Identifier, adapter, registry, and structured-result tests are adopted only
with their corresponding decisions. A narrowed
`test-symbol-ownership.R` should assert one authoritative top-level function
definition after correction consolidation; do not copy its rejected R 4.6,
dependency-removal, exact-five-step, or blanket namespace assumptions.
Provider test families remain companion coverage and are explicitly excluded
from core.

### Vignettes

<!-- split-path status=D path=vignettes/batch_evaluation_metrics.Rmd class=excluded -->
<!-- split-path status=M path=vignettes/proBatch.Rmd class=equivalent -->
<!-- split-path status=M path=vignettes/proBatchFeatures.Rmd class=equivalent -->
<!-- split-path status=D path=vignettes/prone_with_probatch.Rmd class=excluded -->

The two retained core vignettes are byte-identical to the stopped comparator.
Benchmark metrics and PRONE guidance belong to the companion package. Manual
Roxygen/documentation regeneration remains a maintainer task outside this
generated-file-excluding report.

## Companion API disposition

All 26 BEC-only exports are absent from both the stopped core split and the
current destination. They require no core compatibility wrappers because none
belonged to the baseline public API.

<!-- bench-api symbol=correct_with_BERT disposition=companion -->
<!-- bench-api symbol=correct_with_NormAE disposition=companion -->
<!-- bench-api symbol=correct_with_PLSDA_batch disposition=companion -->
<!-- bench-api symbol=correct_with_RUVIII_C disposition=companion -->
<!-- bench-api symbol=correct_with_mComBat disposition=companion -->
<!-- bench-api symbol=correct_with_omicsGMF disposition=companion -->
<!-- bench-api symbol=estimate_omicsGMF_rank disposition=companion -->
<!-- bench-api symbol=imputeMissForest disposition=companion -->
<!-- bench-api symbol=imputeMissForest.ProBatchFeatures disposition=companion -->
<!-- bench-api symbol=imputeMissForest_df disposition=companion -->
<!-- bench-api symbol=imputeMissForest_dm disposition=companion -->
<!-- bench-api symbol=missForestImpute disposition=companion -->
<!-- bench-api symbol=imputePRONE disposition=companion -->
<!-- bench-api symbol=imputePRONE_df disposition=companion -->
<!-- bench-api symbol=imputePRONE_dm disposition=companion -->
<!-- bench-api symbol=impute_and_correct_with_omicsGMF disposition=companion -->
<!-- bench-api symbol=impute_with_omicsGMF disposition=companion -->
<!-- bench-api symbol=impute_with_omicsGMF.ProBatchFeatures disposition=companion -->
<!-- bench-api symbol=calculate_classification_metrics disposition=companion -->
<!-- bench-api symbol=calculate_variance_partition disposition=companion -->
<!-- bench-api symbol=prepare_variance_partition_df disposition=companion -->
<!-- bench-api symbol=plot_variance_partition disposition=companion -->
<!-- bench-api symbol=plot_variance_partition.df disposition=companion -->
<!-- bench-api symbol=plot_intragroup_variation disposition=companion -->
<!-- bench-api symbol=plot_TSNE disposition=companion -->
<!-- bench-api symbol=plot_UMAP disposition=companion -->

Companion dependency ownership includes BERT/BiocParallel integration,
PLSDAbatch/mixOmics, NormAE/Python, RUVIIIC, PRONE, missForest,
omicsGMF/sgdGMF/SingleCellExperiment, variancePartition, Rtsne, umap, plotly,
and benchmark orchestration. Core's independent uses of `sva`,
`SummarizedExperiment`, and `S4Vectors` remain valid.

## Verification and ordering

The required lineage and undeclared-test-dependency fixes come first. Metadata
and release notes follow explicit release/dependency decisions. Correction
deduplication and its narrowed ownership invariant are a separate decision.
Governance and citation recommendations are independent. Identifier policy
precedes registry, structured-result, and matrix-adapter decisions; the
remaining-change plan records the exact file order and verification for each
actionable path.
