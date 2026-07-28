# BEC E. coli Core Migration — Open Questions

## Purpose

This document records issues found during the BEC E. coli core migration that
**no workflow in the ordered migration resolves**. Each entry needs a manual
maintainer decision after the migration finishes.

Workflows `002` through `032` were scanned for coverage of every entry below;
none of them touch these artifacts or dependency names.

## How to use this document

- Add an entry when a review finds a defect or inconsistency that falls outside
  the scope of the remaining ordered workflows.
- Keep one `<!-- open-question id=<kebab-id> status=<open|resolved> area=<area> -->`
  marker per entry so entries stay countable.
- Remove an entry once the maintainer has resolved it, so every entry left in
  this document still needs a decision. The resolving commit is the record.

## Open questions

### Undeclared `reshape2` used by a focused test

<!-- open-question id=reshape2-undeclared status=open area=dependencies -->

- Found: `002_foundation_core_api` source review (2026-07-28).
- Artifacts: `./DESCRIPTION`, `./tests/testthat/test-correct_batch_effect.R`.
- Evidence: the foundation port removed `reshape2` from `Imports`, but the same
  port added a call to `reshape2::dcast()` in the correction test. `reshape2`
  now appears in neither `Imports` nor `Suggests`.
- Effect: `R CMD check` reports an undeclared package used in tests.
- Coverage: no later workflow reads or edits
  `tests/testthat/test-correct_batch_effect.R` for dependency purposes, and no
  workflow mentions `reshape2`.
- Options: (a) add `reshape2` to `Suggests`; (b) rewrite the test to use an
  already-declared package such as `tidyr` or `data.table`; (c) build the matrix
  with the package's own `long_to_matrix()`.
- Decision: _pending_.

### `plyr` declared in `DESCRIPTION` but absent from `pixi.toml`

<!-- open-question id=plyr-missing-from-pixi status=open area=environment -->

- Found: `002_foundation_core_api` source review (2026-07-28).
- Artifacts: `./DESCRIPTION`, `./pixi.toml`.
- Evidence: the port added `plyr` to `Imports` for `plyr::arrange()` in
  `R/plot_split_violinplot.R`, but `pixi.toml` has no `r-plyr` entry. The
  package currently resolves only transitively through `r-reshape2`.
- Effect: the development environment silently depends on a transitive package;
  removing `r-reshape2` from `pixi.toml` would break the build.
- Coverage: no workflow references `pixi.toml`, and no specification governs it.
- Options: (a) add an explicit `r-plyr` pin to `pixi.toml`; (b) replace the
  single `plyr::arrange()` call with `dplyr::arrange()` and drop `plyr` from
  `Imports`.
- Decision: _pending_.

### `pixi.toml` still pins `r-reshape2` after the package dependency was dropped

<!-- open-question id=pixi-stale-reshape2 status=open area=environment -->

- Found: `002_foundation_core_api` source review (2026-07-28).
- Artifacts: `./pixi.toml`, `./DESCRIPTION`.
- Evidence: `reshape2` is no longer a package dependency, but `pixi.toml` still
  pins `r-reshape2 = ">=1.4.5,<2"`.
- Effect: the environment and the package manifest disagree. Resolution depends
  on the outcome of `reshape2-undeclared`.
- Coverage: no workflow references `pixi.toml`.
- Decision: _pending; resolve together with `reshape2-undeclared`._

### `pixi.toml` is outside the specification and Depmesh model

<!-- open-question id=pixi-ungoverned status=open area=governance -->

- Found: `002_foundation_core_api` source review (2026-07-28).
- Artifacts: `./pixi.toml`, `./depmesh.toml`, `./specs/behavior/files_relations.md`.
- Evidence: `depmesh.toml` governs `DESCRIPTION`-adjacent artifacts only through
  R source, test, vignette, workflow, `.Rbuildignore`, `depmesh.toml`, and
  `donna.toml` rules. `pixi.toml` and `DESCRIPTION` have no relation, so a
  dependency change in one is never surfaced when the other is edited.
- Effect: dependency drift between the package manifest and the development
  environment is invisible to Depmesh queries.
- Options: (a) add a Depmesh relation pairing `DESCRIPTION` with `pixi.toml`;
  (b) accept `pixi.toml` as maintainer-owned and document that explicitly.
- Decision: _pending._

### Focused tests without a Depmesh `tests` relation

<!-- open-question id=depmesh-unmapped-tests status=open area=governance -->

- Found: `002_foundation_core_api` source review (2026-07-28).
- Artifacts: `./tests/testthat/test-ProBatchFeatures_links.R`,
  `./tests/testthat/test-explained_variance_plots.R`,
  `./R/batch_correction_helpers.R`, `./R/pbf_input_helpers.R`,
  `./tests/testthat/helper-example-data.R`, `./depmesh.toml`.
- Evidence: `depmesh -p llm deps` returns only `governed_by` for these files.
  `test-ProBatchFeatures_links.R` covers `R/ProBatchFeatures.R` but its filename
  does not match the conventional `test-{module}.R` mapping.
  `test-explained_variance_plots.R` has no matching `R/explained_variance_plots.R`
  because the retained core PVCA behavior lives in `R/proteome_wide_diagnostics.R`.
  `R/batch_correction_helpers.R` and `R/pbf_input_helpers.R` have no `tested_by`
  edge. `helper-example-data.R` is a testthat helper and is not expected to map.
- Effect: dependency discovery under-reports test coverage for these sources.
- Coverage: no later workflow edits `depmesh.toml`.
- Options: (a) add explicit `one_of`/`list` rules to `depmesh.toml` as already
  done for `test-correct_batch_effect.R` and `test-batch_effect_steps.R`;
  (b) rename the test files to match the conventional mapping.
- Decision: _pending._

### Breaking identifier contract has no `NEWS` entry or release-version decision

<!-- open-question id=breaking-contract-news-and-version status=open area=release -->

- Found: `002_foundation_core_api` split review (2026-07-28).
- Artifacts: `./NEWS`, `./DESCRIPTION`.
- Evidence: commit `1adb163d10785e3f775411ec6010b2fb9cd458f8` made duplicate,
  `NA`, and empty sample and feature identifiers a hard error in
  `ProBatchFeatures()`, `long_to_matrix()`, `matrix_to_long()`,
  `.align_sample_annotation()`, and `check_sample_consistency()`. Such
  identifiers were previously repaired silently with `make.unique()` or reported
  only as a warning. The `NEWS` section for `v2.1.0` still describes only the
  design-diagnostics additions and the `center_feature_batch()` deprecation, and
  `DESCRIPTION` keeps `Version: 2.1.0`.
- Effect: an upgrading user hits an error on input that previously loaded, with
  no release note and no version signal that a public contract changed. Later
  children in this migration are expected to add further breaking changes to the
  same unreleased version.
- Also missing after `006`, `009`, `010`, and `011` (2026-07-28): the same
  `v2.1.0` section documents none of the migrated public APIs, including
  `plot_grouped_NA_heatmap()`, the `force_binarization` heatmap threshold, the
  grouped-density `sample_annotation`/`color_by`/`col_vector`/`color_scheme`
  parameters of `plot_NA_density()`, and `pb_subset_samples()`.
- Coverage: no workflow from `003` through `032` reads or edits `NEWS`; only the
  manifest and `002` mention it. `022_72d11d1f7f92_version_variancepartition` is
  the only later workflow that edits `DESCRIPTION`, and it is limited to
  confirming that `2.1.0` already matches the source hunk and is forbidden from
  altering the version without explicit release authority.
  `032_residual_split_review` writes only
  `bec-ecoli-core-residual-review.md` and
  `bec-ecoli-core-remaining-change-plan.md`.
- Options: (a) add a `NEWS` breaking-change block once the migration completes
  and keep `2.1.0`; (b) add the `NEWS` block and bump to a version that signals
  the incompatible contract; (c) relax the new validation to a warning for the
  cases that were previously repaired silently.
- Decision: _pending; requires release authority._

### Shadowed duplicate correction definitions in `R/correct_batch_effects.R`

<!-- open-question id=correction-duplicate-definitions status=open area=api -->

- Found: `003_42544a21f10c_pvca_impl_dedup` split review (2026-07-28), while
  evaluating the comparator's global single-definition test.
- Artifacts: `./R/correct_batch_effects.R`.
- Evidence: `correct_batch_effects_df`, `correct_batch_effects_dm`, and
  `correct_with_removeBatchEffect_dm` are each defined twice at top level in the
  same file, at lines `1169`/`1652`, `1288`/`1689`, and `1367`/`1740`. The
  foundation port collated the historical `R/correct_batch_effects_old.R` bodies
  after the earlier definitions, so the later copies silently win by evaluation
  order. These are the only duplicated top-level function symbols in `R/`; the
  PVCA family is now consolidated.
- Effect: the earlier bodies are unreachable dead code, and the effective
  behavior of three public correction entry points depends on position in a
  1755-line file rather than on an explicit decision. A later migration child
  that edits one of these functions may edit the shadowed copy and observe no
  change. This is the same defect class that source commit `42544a21f10c` fixed
  for PVCA. `R CMD check` does not report it.
- Coverage: the comparator's
  `tests/testthat/test-symbol-ownership.R` asserts that no top-level symbol is
  defined twice, but it was not adopted in `003` because the assertion is
  repository-wide and outside that child's PVCA scope. Of the remaining
  workflows, only `013_fdafc6c8fe13_centralize_all_na_handling` and
  `016_4285c42f3167_clarify_removebatcheffect_missing_policy` touch these
  functions, and both are scoped to missing-value semantics rather than to
  removing the duplicates. `032_residual_split_review` writes reports only.
- Options: (a) delete the shadowed earlier definitions once the intended
  authority is confirmed, and adopt the comparator's repository-wide
  single-definition test; (b) keep the duplicates until `013` and `016` run and
  resolve them there, accepting the risk of editing the shadowed copy;
  (c) resolve after the migration as a separate cleanup.
- Decision: _pending._

### Bayesian PCA and `pcaMethods` skipped by an explicit core-policy decision

<!-- open-question id=bpca-pcamethods-skipped status=open area=api -->

- Found: `011_6d3c15e41905_subset_and_bpca` source review (2026-07-28).
- Artifacts: `./R/proteome_wide_diagnostics.R`, `./DESCRIPTION`,
  `./tests/testthat/test-proteome_wide_diagnostics.R`.
- Evidence: source commit `6d3c15e419058d264b61ac43c27e95e0d635ca01` adds
  `.pb_compute_pca_embedding()` and `.pb_format_pc_axis_label()`, the public
  `pca_method`, `bpca_nPcs`, `bpca_center`, and `bpca_scale` parameters of
  `plot_PCA.default()`, and optional `pcaMethods` in `Suggests`. Slot `011`
  deliberately ported none of it: omitting `fill_the_missing` preserves missing
  values for BPCA while explicitly supplying its existing default `-1` pre-fills
  them, so nominally equivalent calls differ only by argument-supply state; the
  proposed API hides imputation and the fitted model behind four
  backend-specific parameters; and the only source test is conditional on an
  optional dependency that the project Pixi environment does not provide and
  asserts the plot surface rather than the fit.
- Effect: `plot_PCA()` stays SVD-only and still requires complete matrices, so
  missing-aware PCA remains unavailable in lean core. The skip is a reviewed
  policy choice, not an oversight, but it was recorded only in gitignored
  session notes.
- Coverage: no workflow from `012` through `031` mentions `plot_PCA`, BPCA, or
  `pcaMethods`, and `032_residual_split_review` audits the stopped split, which
  is also SVD-only.
- Options: (a) keep core SVD-only; (b) design an explicit missing-policy-aware
  embedding API with fit-level tests and adopt `pcaMethods` in `Suggests`;
  (c) port the source implementation unchanged and accept the
  argument-supply-dependent missing-value semantics.
- Decision: _pending._

## Environment limitations (informational)

These are sandbox properties, not project defects. They are recorded so a future
review does not misread them as regressions.

- `preprocessCore::normalize.quantiles()` cannot create a native worker thread in
  the agent sandbox (`pthread_create()` code 22). Quantile normalization tests
  must be run by the maintainer outside the sandbox.
- `BiocParallel` socket discovery is unavailable in the agent sandbox.
  Registering `BiocParallel::SerialParam()` makes the ComBat-focused tests run.
