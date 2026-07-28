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
- Set `status=resolved` and record the resolution instead of deleting an entry.

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

### Roxygen documentation and `NAMESPACE` regeneration owed

<!-- open-question id=roxygen-regeneration status=open area=documentation -->

- Found: `002_foundation_core_api` source review (2026-07-28); blocking symptom
  identified during the split review (2026-07-28).
- Artifacts: `./R/design_diagnostics.R`, `./R/metadata_diagnostics.R`,
  `./R/proBatch.R`, `./R/correlation-based_diagnostics.R`, and the other ported
  sources.
- Evidence: the migration adds exported functions and new import declarations in
  Roxygen comments, but `man/**` and `NAMESPACE` are maintainer-owned generated
  output that agents must never write. Commit
  `f4d689651c2c3c0e09b2d724db9ca614afb42049` added both the
  `@importFrom tidyr complete nest unnest pivot_longer` declaration in
  `R/proBatch.R` and the unqualified `pivot_longer()` calls in
  `R/correlation-based_diagnostics.R`, so the generated namespace predates the
  import it needs.
- Effect: generated documentation and the namespace lag behind the sources until
  the maintainer regenerates them. This already blocks the test suite:
  `tests/testthat/test-correlation_based_diagnostics.R` reports nine errors of
  the form `could not find function "pivot_longer"`. The failure is stale
  generated output rather than a source defect, and a later review must not read
  it as a regression. Confirmed against an unmodified `HEAD` worktree on
  2026-07-28: the same nine errors occur without any working-tree change.
- Also owed after `003_42544a21f10c_pvca_impl_dedup` (2026-07-28):
  `R/proteome_wide_diagnostics.R` documents the newly activated `stacked_bar`,
  `stacked_plot_title`, `sort_stacked`, `category_order`, `path_to_save_results`,
  and `add_values` parameters, so the PVCA manual pages are stale as well.
- Also owed after the `003` split stage (2026-07-28): the same file now carries
  an explicit `@method plot_PVCA df` tag and a `@details` block beside the
  retained `@export` and `@rawNamespace export(plot_PVCA.df)` tags. The S3
  registration and usage section of the `plot_PVCA.df` manual page therefore
  need regeneration before the generated output can be trusted to match the
  source.
- Coverage: `017_4540aca9182c_generated_missing_docs` is the explicit
  generated-only exception in the family, but regeneration itself stays manual.
- Decision: _pending; regenerate with devtools once the migration completes._

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

## Environment limitations (informational)

These are sandbox properties, not project defects. They are recorded so a future
review does not misread them as regressions.

- `preprocessCore::normalize.quantiles()` cannot create a native worker thread in
  the agent sandbox (`pthread_create()` code 22). Quantile normalization tests
  must be run by the maintainer outside the sandbox.
- `BiocParallel` socket discovery is unavailable in the agent sandbox.
  Registering `BiocParallel::SerialParam()` makes the ComBat-focused tests run.
