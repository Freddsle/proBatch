# BEC E. coli Core Migration History

<!-- bec-ecoli-core-migration-history status=complete slots=32 entries=31 executable-children=24 references=7 residual-paths=99 commits=42 -->

## Status and purpose

The BEC E. coli core migration is complete. This is a non-executable historical record that replaces the retired parent, foundation, per-commit, and residual Donna workflows; their ordered manifest; their three reports; and their temporary progress evidence. It records what was migrated, what was deliberately excluded, which destination commits were accepted, what verification ran, and which post-migration work was transferred to the successor standalone-core baseline workflow.

The migration catalog was marked complete by destination commit `ff5288feb50b47573336e98c558a4a931d2403ac`. The final maintainer-owned Roxygen regeneration was destination commit `968914f264fee562c43a6d321ee36869bdbce639`.

## Immutable provenance

- Source repository: `/home/yuliya/repos/other/proBatch`.
- Source baseline: `ba6ee246eace090e71baa7aba302ca64e76ddb32`.
- Synchronizing merge and foundation snapshot: `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`.
- Source tip: `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.
- Stopped split repository: `/home/yuliya/repos/other/proBatch-core-split`.
- Stopped split preflight: `49cee7cc978fbb149c262a5a783face32dd1d135`.
- Stopped split implementation comparator: `29a7478dc7deea846a2c1ff1abd25a881e6f87db`.
- Stopped split provenance note: `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md`.
- Selected source history: the foundation snapshot followed by the exact 29-commit result of `git rev-list --reverse --topo-order 5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab..e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.
- Ordering: all 29 selected commits are single-parent and linear. The first parent is the foundation merge, each later parent is the preceding selected commit, and the sequence ends at the pinned source tip.
- Synchronizing merge parents: `5b8579bc087396c94de2ed0da21fba34eea81106` and baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32`.
- Synchronizing merge date and message: `2026-03-12T16:12:32+01:00`, `Merge remote-tracking branch 'origin/main' into BEC_ecoli_data`.
- Generated `man/**` and generated `NAMESPACE` were excluded from all agent comparisons. Documentation evidence was taken only from editable Roxygen2 sources under `R/`.

The stopped split was historical comparison evidence, not an accepted patch. Its provenance note says the implementation was rejected and must be reassessed from the baseline. The mode-only preflight changed `DESCRIPTION` from mode `100755` to `100644`. The non-generated baseline-to-split comparison contains 47 paths, comprising 19 additions and 28 modifications, with 9,489 insertions and 677 deletions. It mixed migration with a provider-aware registry, stricter identifier architecture, a matrix adapter, structured step results, dependency removal, and broad lineage changes.

## Ownership boundary at migration closure

The ownership analysis assigned 107 pre-split exports to core. The foundation contained 104; `plot_grouped_NA_heatmap`, `pb_subset_samples`, and `plot_NA_intensity` arrived later in the selected history.

The foundation retained these API groups:

- Container and pipeline: `ProBatchFeatures`, `ProBatchFeatures_from_long`, `as_ProBatchFeatures`, `pb_add_level`, `pb_aggregate_level`, `pb_as_long`, `pb_as_wide`, `pb_assay_matrix`, `pb_current_assay`, `pb_eval`, `pb_pipeline_name`, `get_chain`, `get_operation_log`, `pb_register_step`, `pb_list_steps`, `pb_has_step`, and `pb_transform`.
- Conversion: `long_to_matrix`, `matrix_to_long`, and `create_peptide_annotation`.
- Missing values: `handle_missing_values`, `pb_filterNA`, `pb_groupfilterNA`, `pb_infIsNA`, `pb_nNA`, and `pb_zeroIsNA`.
- Metadata and design: `check_sample_consistency`, `convert_annotation_classes`, `date_to_sample_order`, `dates_to_posix`, `define_sample_order`, `detect_nested_batches`, `detect_outlier_samples`, `filter_metadata_columns`, `find_duplicated_columns`, `guess_factor_columns_if_needed`, `handle_factor_numeric_overlap`, `merge_rare_levels`, `metadata_column_summary`, `subbatch_detection`, `summarize_design`, `validate_batch_design`, and `warn_unmapped_columns`.
- Transformation and normalization: `adjust_batch_trend_df`, `adjust_batch_trend_dm`, `center_feature_batch`, `center_feature_batch_means_df`, `center_feature_batch_means_dm`, `center_feature_batch_medians_df`, `center_feature_batch_medians_dm`, `fit_nonlinear`, `log_transform_df`, `log_transform_dm`, `normalize_data_df`, `normalize_data_dm`, `normalize_sample_medians_df`, `normalize_sample_medians_dm`, `quantile_normalize_df`, `quantile_normalize_dm`, `unlog_df`, and `unlog_dm`.
- Batch correction: `correct_batch_effects`, `correct_batch_effects_df`, `correct_batch_effects_dm`, `correct_with_ComBat`, `correct_with_ComBat_df`, `correct_with_ComBat_dm`, `correct_with_removeBatchEffect`, `correct_with_removeBatchEffect_df`, and `correct_with_removeBatchEffect_dm`.
- Numerical diagnostics and PVCA: `calculate_feature_CV`, `calculate_peptide_corr_distr`, `calculate_sample_corr_distr`, `calculate_PVCA`, and `prepare_PVCA_df`.
- Visualization: `plot_CV_distr`, `plot_CV_distr.df`, `plot_NA_density`, `plot_NA_frequency`, `plot_NA_heatmap`, `plot_PCA`, `plot_PVCA`, `plot_PVCA.df`, `plot_PVCA_stacked_from_saved`, `plot_boxplot`, `plot_corr_matrix`, `plot_heatmap_diagnostic`, `plot_heatmap_generic`, `plot_hierarchical_clustering`, `plot_iRT`, `plot_peptide_corr_distribution`, `plot_peptide_corr_distribution.corrDF`, `plot_peptides_of_one_protein`, `plot_protein_corrplot`, `plot_sample_corr_distribution`, `plot_sample_corr_distribution.corrDF`, `plot_sample_corr_heatmap`, `plot_sample_mean`, `plot_single_feature`, `plot_spike_in`, `plot_split_violin_with_boxplot`, `plot_with_fitting_curve`, `generate_colors_for_numeric`, and `sample_annotation_to_colors`.

Provider implementations, provider registrations, benchmark orchestration, optional-provider discovery, Python integration, and companion-only dependencies were excluded even when they shared a source file or commit with a retained core hunk.

### UMAP and t-SNE supersession

At migration closure, `plot_TSNE` and `plot_UMAP` were classified as companion-owned benchmark plots and deliberately excluded. On 2026-07-29 the maintainer superseded that forward-looking ownership decision: both plots are now required in the standalone core baseline so that `proBatchBench` can import the baseline package. This does not retroactively change the completed migration ledger; implementation and dependency choices belong to the successor standalone-core workflow.

The same standalone/import goal requires a fresh decision on the smallest provider-neutral Core extension surface actually used by proBatchBench. That review does not accept the stopped split's registry, matrix-adapter, or structured-result architecture wholesale; any retained capability must be independently specified and tested by the successor workflow.

## Ordered completion ledger

Slot `001` was the parent coordinator. It verified the pinned inputs, created and governed the ordered child suite, recorded all seven reference-only positions, drove the maintainer review gates, ran the residual audit, and closed the catalog in `ff5288feb50b47573336e98c558a4a931d2403ac`.

The outcome totals are:

- `source-and-split`: 7 entries (`002`, `003`, `011`, `013`, `016`, `021`, and `027`);
- `source-only`: 10 entries (`004`, `006`, `007`, `009`, `010`, `012`, `015`, `019`, `020`, and `031`);
- `split-only`: 1 entry (`028`);
- `no-change`: 5 entries (`008`, `014`, `017`, `022`, and `023`);
- `reference-only`: 7 entries (`005`, `018`, `024`, `025`, `026`, `029`, and `030`);
- `reports`: 1 entry (`032`).

Together these are 31 ordered entries: 23 executable foundation/commit entries, seven deterministic reference-only entries, and one residual review. There were 24 executable child workflows and 25 numbered workflow artifacts when the retired parent is included.

| Slot | Immutable source | Outcome | Destination source commit | Destination split/report commit | Accepted result |
|---|---|---|---|---|---|
| 001 | Migration coordination | coordinator | none | `ff5288feb50b47573336e98c558a4a931d2403ac` | Verified the immutable migration scope, coordinated all ordered child/reference outcomes and review gates, accepted the residual audit, and marked the retired workflow family complete. |
| 002 | Foundation `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab` | source-and-split | `f4d689651c2c3c0e09b2d724db9ca614afb42049` | `1adb163d10785e3f775411ec6010b2fb9cd458f8` | Ported the core API and tests, excluded providers, then adopted strict rejection of duplicate, missing, and empty identifiers. |
| 003 | `42544a21f10ca6960d3e4c44d2833f764054d721` | source-and-split | `8a93e6451194abc95444b3231852f39459d02b94` | `8cc4cb675a3239c74641e36423e473be88241409` | Consolidated PVCA implementations and declared `plot_PVCA.df` as an S3 method. |
| 004 | `db128668458a58ad31f66be5b6e39e2fedadbbe1` | source-only | `a29b6db34f587f0b2304cdd5edc249c8672ed941` | none | Retained only the core QFeatures link guard and focused duplicate/missing feature-axis tests; excluded PRONE integration. |
| 005 | `65f70a46c4cf44e2717744aeafbf8acbe83b0378` | reference-only | none | none | Recorded PRONE condition/batch metadata, backend compatibility, provider tests, and vignette as companion-only. |
| 006 | `4e4e9811503cc7da1e35a8113a2fb8383fc5007b` | source-only | `6d59081b4b4b52501ff8635b3d6e7285403d2b55` | none | Added grouped missingness heatmaps, annotation alignment, validation, and multi-assay coverage. |
| 007 | `12d4597e9e86ea2daddc77c4d236e692bbfa8fa9` | source-only | `6bbbee20eee64efe511e83db9014c28bf3bda0ba` | none | Qualified split-violin rendering calls and kept `variancePartition` outside core. |
| 008 | `5747090a1de1cdde780c6db7a84f988aedb9e8af` | no-change | none | none | Test-only binarization changes intentionally awaited slot 009; later focused coverage passed. |
| 009 | `7feb6d5cea9ec92d378453ff65d672f231fa44b8` | source-only | `2029fe0e65d198594226c72e1c9c89d1c8959e92` | none | Added grouped NA-density plots and grouped-heatmap binarization with focused tests. |
| 010 | `0dbce8129444c9079a560fb9bcd60b328b1f054c` | source-only | `4e91d9ce39d7929bca3c6d1b96cfa432fd7cbd5c` | none | Added named and palette grouped-density color schemes and tests. |
| 011 | `6d3c15e419058d264b61ac43c27e95e0d635ca01` | source-and-split | `a590480a5137cd6df40421748ca17d565cddfd18` | `2d3026e5ccb7d4122410fde0a5c2334aea58577c` | Added sample subsetting, correct four-index dispatch, and per-assay vector/title handling; deliberately kept PCA SVD-only and recorded the BPCA decision. |
| 012 | `7db05faed18f75572d92553ae668288c326bda2a` | source-only | `9e3ef760e24ef10e21d9e7510263d4fb074a19d0` | none | Accepted symmetric integer-like factor/numeric `colData` compatibility while rejecting conversion-created NA matches and excluding omicsGMF/PRONE. |
| 013 | `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5` | source-and-split | `d66b2bcdab3c37d140b532a9e435805fb098643e` | `53efd441b5f068b9c0c82d0b337d3f48d758a067` | Hardened wide correction inputs/defaults and preserved explicit `removeBatchEffect` forwarding; provider-only all-NA helpers stayed excluded. |
| 014 | `77a3365f0e13e9abca40a14e1898a351a1cdbc7f` | no-change | none | none | The provider-only all-NA helpers had no destination target, so their Roxygen demotion required no change. |
| 015 | `4b117a03f0c5295f2253db26e13b945d7e3027b5` | source-only | `52704e234ee2ad85e33f24c28886a8c1eff7605d` | none | Made failing-group masking default-on and added `plot_NA_intensity` for matrices, `SummarizedExperiment`, and single-/multi-assay `ProBatchFeatures`. |
| 016 | `4285c42f31670d2f750dc8eb8c7ff1d0134a342d` | source-and-split | `5b034afddef2c134c071a74627ca8d32dff09686` | `21cab206a68ffdcec75cb51319f662829c62030f` | Clarified legacy removeBatchEffect behavior, then standardized correction missing policies as `error`, `keep`, `drop_features`, and `fill` with `fill_value`. |
| 017 | `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57` | no-change | none | none | Confirmed a generated-only documentation commit had no editable target without inspecting generated documentation. |
| 018 | `b8a262b4256966d60e1e8452ebde7a1bf471b4af` | reference-only | none | none | Recorded omicsGMF batch-mean correction semantics and provider tests as companion-only. |
| 019 | `d95dd736cb27d68fb2e21b20ab97d0f69826a663` | source-only | `9fbe2bc8e4174077a330c06babb9b0b2ec4fb479` | none | Made zero/near-zero/non-finite CV undefined, restored date locale, made tie order deterministic, masked LOESS extrapolation, and excluded mComBat. |
| 020 | `6989542e4ece84f9fe210f3bf06b8254684e425e` | source-only | `f69bdc029be8e2ff2af4cd38d0d70805dd241bda` | none | Hardened pairwise correlations, nonlinear copy safety, detailed sample-consistency warnings, and correction survivor ordering; excluded PRONE. |
| 021 | `428ad74b73fbc067de09d0a71c93395ae85eb51f` | source-and-split | `2a60cbef23b45a6d3edd52a3937f9dd539b6c67c` | `86fac7c5bb538563135f201feadf1f3e7c150073` | Localized LOESS warning capture, documented core returns, and made NA-intensity examples self-contained; provider and embedding hunks stayed excluded. |
| 022 | `72d11d1f7f92a10cb9a5244e68f9c19e830197b9` | no-change | none | none | Kept destination version `2.1.0` and omitted companion-owned `variancePartition`; neither source nor stopped split was release authority. |
| 023 | `8088de1701e0908cca25b978105f8d6c7bfccc20` | no-change | none | none | Kept existing Docker ignore policy and rejected speculative `Docker_README.md`/`Rplots.pdf` ignore changes. |
| 024 | `4a06d99949114b2804b8d34492a288872fb611ed` | reference-only | none | none | Recorded missing-policy forwarding through BERT, PLSDA, omicsGMF, PRONE, and missForest adapters as companion-only. |
| 025 | `6601232c69b44c507f2e9f63d256836e655f7973` | reference-only | none | none | Recorded omicsGMF design-matrix column-name restoration as companion-only. |
| 026 | `96d38eb7449e4d38c0d6a3fffd66e17f145f669f` | reference-only | none | none | Recorded omicsGMF latent artifacts on correction/imputation outputs as companion-only comparison evidence. |
| 027 | `b60e2b169afe45934c18eee1c83469e7e5fed33f` | source-and-split | `0ff70e84956c28bde2996c860b17460d6404faf6` | `13c712caaa2510ba5e811db5ec714a6556918f81` | Materialized explicitly named final log assays and adopted strict final-name validation, assay-specific lineage resolution, and lineage-derived pipeline naming. |
| 028 | `5edcc4b89f5cbfba07e1ca1057ab475a56501202` | split-only | none | `b0c9b04a5dc13cebc52db8b30ad414dba1d32be1` | Kept the evolved distinction between link pruning and warned link removal, adding focused in-place filter coverage without obsolete custom warnings. |
| 029 | `0b3c3b55403f2bb7342a37236dd30f3d9b3544e9` | reference-only | none | none | Recorded the PRONE normalization default-input change to `log2` as companion-only. |
| 030 | `20e76c9b9b28c0ec98faa63c3f9382c1347301b9` | reference-only | none | none | Recorded combined omicsGMF imputation/correction, aliases, and tests as companion-only. |
| 031 | `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92` | source-only | `91fb4dde931c1bb6fe5cc5c9d9eea03f6c9167db` | none | Dropped unused batch and factor-covariate levels before ComBat/removeBatchEffect model construction; excluded BERT. |
| 032 | Residual stopped-split review | reports | none | `ea91575e02cb3a6dab3ceeae9e31c3d29e49c3c3` | Classified the immutable 99-path residual, recorded all companion API dispositions, and produced an ordered follow-up plan. |

Every executable foundation/commit entry passed both maintainer review pauses and both post-resume reverification gates. Every reference entry recorded only its deterministic no-commit disposition.

## Accepted behavior after the ordered migration

The resulting core baseline accepted these material changes:

- The core-owned API at the synchronization snapshot was transferred without provider implementations or provider-only dependencies.
- Duplicate, missing, and empty feature/sample identifiers became hard errors across container, conversion, annotation-alignment, and consistency entry points.
- PVCA implementation ownership was consolidated and its data-frame plot registered as an S3 method.
- QFeatures assay links are created only when feature axes permit a valid one-to-one link.
- Grouped missingness heatmaps, binarization, grouped density plots, color schemes, and NA-intensity plots became core plotting behavior.
- Sample subsetting preserves assay structure, supports four-index dispatch, and handles shared or per-assay atomic vectors explicitly.
- Matching integer-like factor and numeric sample metadata are accepted symmetrically without accepting genuine conflicts or conversion-created missing values.
- Wide correction wrappers validate numeric matrices, sample names, annotation membership, defaults, and explicit method forwarding.
- Grouped missing filtering masks failing groups by default.
- Correction APIs use the canonical `error`, `keep`, `drop_features`, and `fill` missing policies with explicit `fill_value`.
- CV returns `NA` for zero, near-zero, and non-finite denominators; date conversion restores locale state; tied sample times are deterministic; LOESS corrections do not extrapolate beyond fitted support.
- Sample correlations use pairwise-complete observations, nonlinear fitting protects caller-owned inputs, warnings identify mismatched samples precisely, and correction restores surviving sample order.
- Explicit final log names materialize transformed assays, final assay names are validated, and pipeline/lineage resolution is assay-specific.
- In-place missing filters have focused coverage for their distinct QFeatures link effects.
- ComBat and removeBatchEffect drop unused model levels before fitting.
- Editable Roxygen sources were regenerated by the maintainer at the end of the migration, making `plot_NA_intensity` and its S3 methods reachable.

## Deliberate no-change, reference, and exclusion decisions

The five no-change entries were not omissions:

- Slot `008` preserved source order while slot `009` supplied the implementation and end-to-end tests.
- Slot `014` targeted provider-only helpers that were not retained.
- Slot `017` contained generated documentation only; editable documentation was already owned by slot `016`.
- Slot `022` did not treat either historical version metadata or `variancePartition` as core release authority.
- Slot `023` did not add speculative ignore rules.

The seven reference-only entries retained immutable evidence for companion changes without creating source or split commits: PRONE metadata (`005`), omicsGMF batch-mean semantics (`018`), provider missing-policy forwarding (`024`), omicsGMF design names (`025`), omicsGMF latent artifacts (`026`), PRONE default input (`029`), and combined omicsGMF imputation/correction (`030`).

The migration deliberately did not copy the stopped split's broad redesign:

- `pb_unregister_steps` and provider-aware registry metadata, collision, availability, and replay semantics;
- `pb_step_result` and structured artifact/materialization semantics;
- `pb_apply_matrix_method` and its coupled matrix/long adapter policy;
- additional positional-rowname prohibitions, mandatory feature identifiers, and duplicate aggregation policy;
- blanket dependency removal, exact registry-step assertions, and repository-wide ownership assumptions;
- unchanged BPCA code with argument-supply-dependent missing-value behavior and backend-specific public parameters.

The provider correction, imputation, benchmark, orchestration, Python, and vendored resource families remained outside core. At migration closure the 26 companion exports were:

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

The final two markers describe the historical migration decision only; the dated UMAP/t-SNE supersession above governs forward planning.

## Residual stopped-split audit

The exhaustive non-generated comparison `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92..29a7478dc7deea846a2c1ff1abd25a881e6f87db` contains exactly 99 paths: 11 added, 38 deleted, and 50 modified. Its exact name-status stream has SHA-256 `28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38`.

The reviewed classifications were:

- `required`: 4 paths;
- `recommended`: 7 paths;
- `decision`: 13 paths;
- `equivalent`: 43 paths;
- `excluded`: 32 paths.

`Equivalent` includes destination behavior stronger than the stopped split. `Excluded` identifies companion-owned implementation, tests, resources, and vignettes. `Required`, `recommended`, and `decision` described post-migration work rather than incomplete execution of the migration itself.

The exact ordered path/class evidence is:

<!-- split-path status=M path=.Rbuildignore class=equivalent -->
<!-- split-path status=M path=.gitignore class=equivalent -->
<!-- split-path status=M path=DESCRIPTION class=decision -->
<!-- split-path status=M path=NEWS class=required -->
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
<!-- split-path status=M path=README.md class=equivalent -->
<!-- split-path status=M path=inst/CITATION class=recommended -->
<!-- split-path status=D path=inst/overrides/PRONE/Imputation.R class=excluded -->
<!-- split-path status=D path=inst/overrides/normae/__main__.py class=excluded -->
<!-- split-path status=D path=inst/scripts/README.md class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_task_helpers.R class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_workflow.R class=excluded -->
<!-- split-path status=D path=inst/scripts/batch_correction_workflow_params.R class=excluded -->
<!-- split-path status=D path=inst/scripts/pb_tasks.yaml class=excluded -->
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
<!-- split-path status=D path=vignettes/batch_evaluation_metrics.Rmd class=excluded -->
<!-- split-path status=M path=vignettes/proBatch.Rmd class=equivalent -->
<!-- split-path status=M path=vignettes/proBatchFeatures.Rmd class=equivalent -->
<!-- split-path status=D path=vignettes/prone_with_probatch.Rmd class=excluded -->

## Follow-ups transferred at closure

The completed migration intentionally did not resolve every defect or release task discovered during review. The successor standalone-core baseline workflow owns the following work; this historical list is evidence, not an executable plan:

| Area | Evidence at closure | Transferred work |
|---|---|---|
| Lineage consistency | `pb_add_level` accepted an existing target without checking requested data/parent lineage, and the operation log permitted a conflicting second non-self parent. | Add narrow conflict validation and focused regression tests without importing the rejected registry architecture. |
| Missing-filter output names | Explicit `final_name` collisions in `pb_filterNA` and `pb_groupfilterNA` could be silently suffixed. | Define and test a consistent explicit-name collision error while preserving justified generated-name behavior. |
| Focused-test dependency | `test-correct_batch_effect.R` used undeclared `reshape2::dcast`. | Replace it with an independent base-R construction and remove stale `r-reshape2` environment support. |
| Runtime dependency | The sole `plyr` use was `plyr::arrange` in the split-violin implementation, while Pixi supplied it only transitively. | Replace it with stable base ordering, remove `plyr`, and synchronize package/environment metadata. |
| Build boundary | `.pixi`, `pixi.toml`, and `pixi.lock` were not excluded by `.Rbuildignore`; `.pixi` was approximately 2.2 GB. | Add anchored build exclusions before any package archive or broader check. |
| Test framework | Tests use `local_mocked_bindings` and `expect_no_warning`, newer than the declared `testthat (>= 3.0.0)`. | Raise the minimum to at least `3.2.0` and synchronize the development environment. |
| Correction ownership | `correct_batch_effects_df`, `correct_batch_effects_dm`, and `correct_with_removeBatchEffect_dm` each had a shadowed earlier definition. | Retain the effective later deprecated forwarding wrappers, delete dead copies, and add a narrow top-level uniqueness invariant. |
| Test simplification | `test-ProBatchFeatures_links.R` duplicated the authoritative container suite. | Merge any unique value, then remove the redundant file rather than govern duplicate coverage. |
| Depmesh coverage | Helper and non-conventional source/test pairs were unmapped. | Map batch-correction helpers to correction tests, PBF input helpers to plot-helper tests, and proteome-wide diagnostics to explained-variance tests. |
| Release notes | `NEWS` did not cover strict identifiers, missing-policy/default changes, group masking, new plots/subsetting, CV/LOESS changes, lineage, or correction hardening. | Document the accepted public and breaking behavior; keep the version decision under release authority. |
| Vignettes | Core vignettes still taught legacy `fill_the_missing`, `"remove"`, an obsolete `-1` default, and deprecated correction wrappers. | Update to canonical missing policies and unified APIs, then build vignettes in the maintainer environment. |
| Dependency cleanup | `corrplot`, `ggfortify`, and `lazyeval` appeared import-only; `gridExtra` was optional/qualified; `devtools` and `gtable` lacked package use. | Verify each candidate narrowly, synchronize Roxygen imports and manifests, and avoid the stopped split's blanket purge. |
| BPCA | The proposed `pcaMethods` surface had argument-supply-dependent missing semantics and weak fit-level coverage. | Keep SVD-only unless a separate, explicit missing-aware embedding design is approved. |
| UMAP and t-SNE | They were historically excluded as benchmark plots. | Per the 2026-07-29 supersession, migrate suitable original proBatch implementations into the standalone core with focused tests and minimal optional dependencies. |
| Citation | The citation used deprecated construction and misspelled `Evan G. Williams`. | Correct and modernize it as independent package maintenance. |
| Roxygen | The former open question about `plot_NA_intensity` registration was stale after `968914f264fee562c43a6d321ee36869bdbce639`. | No migration action remains; future Roxygen edits still require maintainer regeneration. |
| Rejected split architecture | Registry, adapter, structured-result, and expanded identifier policies were coupled redesigns. | Keep them out of the baseline closure workflow unless separately specified and approved. |

## Verification evidence and limitations

- Donna's completion gate accounted for all 31 ordered entries, the 24 executable children, the seven reserved reference records, their destination commits/no-change decisions, the residual report commit, the 99-path checksum, and the 26 historical companion API dispositions.
- All retired Donna workflow artifacts validated and rendered before or during execution, and all maintained source/split review gates were explicitly resumed by the maintainer.
- Focused workflow evidence included: 59 assertions for the slot-008 grouped-heatmap file; 73 for grouped density/binarization; 83 for density colors; 192 for sample subsetting/PCA-related behavior; 156 for `colData` compatibility; 156 plus 55 for centralized all-NA/correction behavior; 162 for group masking and NA-intensity; 166, 609, and 106 for the two missing-policy stages and diagnostics; 20 for CV/date/LOESS safeguards; 10 for alignment/consistency safeguards; 315 plus 59 for documentation/LOESS and NA-intensity; and 113 plus 28 for unused correction levels.
- The NA-intensity color-focused checks retained two known rank-deficient spline warnings. Their cleanup was transferred to the successor workflow.
- Post-closure focused auditing also observed 147 passing container assertions, 24 passing assertions in the redundant link file, 35 explained-variance assertions, 13 split-violin assertions, 59 NA-intensity assertions with the two known warnings, and 113 correction assertions.
- Maintainer regeneration commit `968914f264fee562c43a6d321ee36869bdbce639` made `plot_NA_intensity` exported and registered its `default` and `ProBatchFeatures` S3 methods.
- No agent read, compared, generated, linted, or validated a file under `man/`. No agent-run package-wide check that could inspect `man/` was used.
- `preprocessCore::normalize.quantiles` could not create a native worker thread in the sandbox. Quantile-normalization verification must run outside it.
- BiocParallel socket discovery was unavailable in the sandbox. Default ComBat behavior must be verified in the maintainer environment; focused tests used an explicit serial registration where necessary.

## Destination commit ledger

These are the exact 42 destination commits from the initial migration workflow through final Roxygen regeneration:

1. `532f42f0cd3ece41c95f4ad56af572736661e09d` — chore(workflows): add BEC E. coli core migration workflow
2. `c29730c2ca6c107e31cd9ad190e9fa4cdeb178af` — feat(pixi): add pixi.lock and pixi.toml for environment management
3. `52cdc74f54e110ad9e6e937460eda388babe4fb4` — chore(workflows): add BEC core migration workflow suite
4. `f4d689651c2c3c0e09b2d724db9ca614afb42049` — feat: port core API from BEC synchronization snapshot
5. `4886d0558514a32e2ef550dea56c1f7551e5be6a` — chore(workflows): govern reports and track migration open questions
6. `1adb163d10785e3f775411ec6010b2fb9cd458f8` — feat!: reject duplicate and invalid identifiers across core entry points
7. `cebd254b906dd1d4f2e76ca1a939a398118765a1` — docs(workflows): record migration release and regeneration open questions
8. `8a93e6451194abc95444b3231852f39459d02b94` — feat(pvca): port upstream consolidation behavior
9. `8cc4cb675a3239c74641e36423e473be88241409` — chore: declare plot_PVCA.df method and record split-review findings
10. `a29b6db34f587f0b2304cdd5edc249c8672ed941` — feat(tests): add tests for handling duplicate and missing feature axes in add-assay-with-link
11. `6d59081b4b4b52501ff8635b3d6e7285403d2b55` — feat(missing-data): add grouped missingness heatmaps
12. `9248ce8a5975b2fcaf8dd996987a8363e165c611` — docs(workflows): record owed regeneration for the grouped NA heatmap API
13. `9fac6acecb445e9fff02892991a48fef898de8f8` — build(package): regenerate Roxygen outputs
14. `0cf2cb8443ba4ed6c9c1f72adbcf2470224b3cb0` — docs(workflows): drop the resolved regeneration open question
15. `7a9c9385c94f8312a5595f062ba50adaaf09466a` — chore: add styler to Pixi and format missing-plot test
16. `6bbbee20eee64efe511e83db9014c28bf3bda0ba` — fix(plots): qualify split-violin rendering helpers
17. `4fd90f249b9903addb6f88e8c6643a0263624999` — chore(workflows): add commit message handoffs at review gates
18. `2029fe0e65d198594226c72e1c9c89d1c8959e92` — feat(missing-data): add grouped density and heatmap binarization
19. `4e91d9ce39d7929bca3c6d1b96cfa432fd7cbd5c` — feat(plots): add grouped density color schemes
20. `a590480a5137cd6df40421748ca17d565cddfd18` — feat(pbf): improve sample subsetting and per-assay titles
21. `2d3026e5ccb7d4122410fde0a5c2334aea58577c` — docs(workflows): record skipped Bayesian PCA and missing NEWS entries
22. `745a211c04fd99677ce35b0a940e85be9712757a` — build(package): regenerate Roxygen outputs
23. `8c1d3e30ae09f891ef14816622ac63278b582adc` — style(tests): reformat grouped density row ordering
24. `9e3ef760e24ef10e21d9e7510263d4fb074a19d0` — fix(pbf): accept matching integer-like factor metadata
25. `d66b2bcdab3c37d140b532a9e435805fb098643e` — fix(batch-correction): validate wide inputs and select wrapper defaults
26. `53efd441b5f068b9c0c82d0b337d3f48d758a067` — fix(batch-correction): allow removeBatchEffect in wide wrapper
27. `52704e234ee2ad85e33f24c28886a8c1eff7605d` — feat(missing-data)!: add failing-group masking and NA-intensity plots
28. `5b034afddef2c134c071a74627ca8d32dff09686` — fix(batch-correction): align removeBatchEffect missing-value handling
29. `dc16daeaa3b158886ae98cea2555f38f3ce44e35` — docs(workflows): refresh migration open questions after slots 013-016
30. `21cab206a68ffdcec75cb51319f662829c62030f` — feat(missing-data)!: standardize batch-correction policies
31. `9fbe2bc8e4174077a330c06babb9b0b2ec4fb479` — fix: stabilize CV, date handling, and LOESS boundaries
32. `ccc5343c5d660cc4e28b9bbcd01f28ca1a65caea` — docs(workflows): update migration open questions after slots 016-019
33. `f69bdc029be8e2ff2af4cd38d0d70805dd241bda` — fix: harden sample consistency, correlations, and nonlinear fitting
34. `2a60cbef23b45a6d3edd52a3937f9dd539b6c67c` — refactor(package): localize LOESS warnings and document returns
35. `86fac7c5bb538563135f201feadf1f3e7c150073` — docs(plots): make NA-intensity examples self-contained
36. `0ff70e84956c28bde2996c860b17460d6404faf6` — fix(pbf): materialize explicitly named final log assays
37. `13c712caaa2510ba5e811db5ec714a6556918f81` — fix(pbf)!: validate final assay names and resolve lineage
38. `b0c9b04a5dc13cebc52db8b30ad414dba1d32be1` — test(missing-data): cover in-place filter link semantics
39. `91fb4dde931c1bb6fe5cc5c9d9eea03f6c9167db` — fix(correction): drop unused factor levels before model fitting
40. `ea91575e02cb3a6dab3ceeae9e31c3d29e49c3c3` — docs(workflows): record residual split audit and change plan
41. `ff5288feb50b47573336e98c558a4a931d2403ac` — docs(workflows): mark BEC migration workflows implemented
42. `968914f264fee562c43a6d321ee36869bdbce639` — build(package): regenerate Roxygen outputs
