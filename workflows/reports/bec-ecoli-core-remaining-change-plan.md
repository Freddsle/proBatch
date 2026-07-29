# BEC E. coli Core Remaining Change Plan

## Purpose

This plan contains one entry for every `required`, `recommended`, or `decision`
path in the residual review. Orders are execution order, not Git inventory
order. A `decision` entry requires maintainer approval before its dependent
implementation proceeds.

## 1. Correctness and focused regression

<!-- plan-path path=R/ProBatchFeatures.R class=required destination=@/R/ProBatchFeatures.R order=1 reason=reject-conflicting-results -->

Reject an existing assay target when requested data or parent lineage differs,
and reject a second non-self parent edge for an existing operation-log result.
No dependency change is required; run the focused container and link suites.

<!-- plan-path path=tests/testthat/test-lineage.R class=required destination=@/tests/testthat/test-ProBatchFeatures.R order=2 reason=lineage-conflict-regressions -->

Selectively merge the existing-target data/lineage mismatch and conflicting
parent-edge cases into the current container suite. Do not import already
equivalent or registry-coupled cases wholesale.

<!-- plan-path path=tests/testthat/test-correct_batch_effect.R class=required destination=@/tests/testthat/test-correct_batch_effect.R order=3 reason=remove-undeclared-reshape2 -->

Replace `reshape2::dcast()` with `long_to_matrix()` or another already-declared
conversion, preserving the ComBat regression. Verify the full correction test
and confirm no test-only undeclared dependency remains.

## 2. Dependency and release decisions

<!-- plan-path path=R/plot_split_violinplot.R class=recommended destination=@/R/plot_split_violinplot.R order=4 reason=remove-plyr-call -->

Replace the sole `plyr::arrange()` call with stable base ordering while
retaining explicit ggplot2 qualification and current plot behavior. Run the
split-violin focused test.

<!-- plan-path path=DESCRIPTION class=decision destination=@/DESCRIPTION order=5 reason=resolve-release-dependencies -->

Decide the post-migration version signal and optional dependency policy. Keep
the supported R version authoritative, remove `plyr` only after order 4, and
do not copy the rejected split's R 4.6, BiocParallel default, or `1.99.5`
metadata. Parse as DCF and reconcile the development environment separately.

<!-- plan-path path=R/proBatch.R class=decision destination=@/R/proBatch.R order=6 reason=synchronize-roxygen-imports -->

If optional dependencies are removed, update only their Roxygen import sources
and qualify or replace every remaining use. Leave generated documentation and
`NAMESPACE` regeneration to the maintainer.

<!-- plan-path path=NEWS class=required destination=@/NEWS order=7 reason=document-breaking-migration -->

After the API and release decisions settle, document every accepted
user-visible and breaking migration change under the chosen release heading.
Review against the confirmed migration commits and public Roxygen sources.

## 3. Correction ownership

<!-- plan-path path=R/correct_batch_effects.R class=decision destination=@/R/correct_batch_effects.R order=8 reason=choose-authoritative-definitions -->

Confirm that the effective later definitions of `correct_batch_effects_df`,
`correct_batch_effects_dm`, and `correct_with_removeBatchEffect_dm` remain
authoritative, then remove only their shadowed copies. Re-run correction and
registered batch-step suites.

<!-- plan-path path=tests/testthat/helper-source-root.R class=recommended destination=@/tests/testthat/helper-source-root.R order=9 reason=centralize-source-root -->

Add one safe test source-root helper instead of duplicating repository lookup.
It has no package runtime dependency and supports the ownership regression.

<!-- plan-path path=tests/testthat/test-symbol-ownership.R class=decision destination=@/tests/testthat/test-symbol-ownership.R order=10 reason=add-narrow-uniqueness-guard -->

If order 8 is approved, add only the top-level single-definition invariant and
the accepted provider boundary. Reject the comparator's exact R version,
dependency set, five-step registry, and blanket namespace assertions.

## 4. Dependency-discovery coverage

<!-- plan-path path=R/batch_correction_helpers.R class=recommended destination=@/depmesh.toml order=11 reason=map-helper-coverage -->

Add the focused correction tests as `tested_by` coverage without inlining or
deleting this retained core helper; verify both Depmesh directions.

<!-- plan-path path=R/pbf_input_helpers.R class=recommended destination=@/depmesh.toml order=12 reason=map-input-helper-coverage -->

Add the current container/input-helper coverage relation and its reverse,
preserving the accepted helper implementation.

<!-- plan-path path=tests/testthat/test-ProBatchFeatures_links.R class=recommended destination=@/depmesh.toml order=13 reason=map-link-tests -->

Map this retained link suite to `R/ProBatchFeatures.R`; do not follow the
stopped split's wholesale test deletion.

<!-- plan-path path=tests/testthat/test-explained_variance_plots.R class=recommended destination=@/depmesh.toml order=14 reason=map-pvca-tests -->

Map retained core PVCA coverage to `R/proteome_wide_diagnostics.R`; companion
variance-partition behavior remains excluded.

## 5. Independent citation maintenance

<!-- plan-path path=inst/CITATION class=recommended destination=@/inst/CITATION order=15 reason=modernize-citation -->

Use `bibentry(bibtype = "Article")` and correct `Evan G. Williams`. Parse the
citation file and confirm the DOI and author list are otherwise unchanged.

## 6. Identifier and extension API decisions

<!-- plan-path path=R/identifiers.R class=decision destination=@/R/identifiers.R order=16 reason=decide-identifier-policy -->

Decide positional annotation row names, mandatory feature identifiers, and
intentional duplicate-key aggregation together. Preserve the already accepted
duplicate and missing identifier errors.

<!-- plan-path path=tests/testthat/test-identifiers.R class=decision destination=@/tests/testthat/test-identifiers.R order=17 reason=cover-identifier-policy -->

Add only cases matching the approved identifier policy, then run conversion,
container, correction-alignment, and consistency suites.

<!-- plan-path path=R/registry.R class=decision destination=@/R/registry.R order=18 reason=decide-provider-registry -->

Decide provider metadata, aliases, collision/overwrite rules, availability
checks, replay semantics, and `pb_unregister_steps()`. Preserve current core
steps, including `loessLimmaRBE`, and require no companion package at core
load time.

<!-- plan-path path=R/zzz_helpers.R class=decision destination=@/R/zzz_helpers.R order=19 reason=align-registry-lifecycle -->

Change startup registration only if order 18 is approved. Keep loading silent,
core-only, and free of package installation or source-time provider discovery.

<!-- plan-path path=tests/testthat/test-registry.R class=decision destination=@/tests/testthat/test-registry.R order=20 reason=verify-provider-registry -->

Adopt collision, availability, replay, and ownership tests only with the
registry contract. Update the exact core-step expectation for all accepted
registrations.

<!-- plan-path path=R/step_result.R class=decision destination=@/R/step_result.R order=21 reason=decide-structured-results -->

Decide artifact names, validation, replacement/merge behavior, metadata
storage, and virtual-step behavior before exposing `pb_step_result()`.

<!-- plan-path path=R/matrix_adapter.R class=decision destination=@/R/matrix_adapter.R order=22 reason=decide-matrix-adapter -->

After identifier, registry, and structured-result decisions, decide
`pb_apply_matrix_method()` semantics for long-row order, ordered subsets,
missing policies, annotation supply, and `keep_all`.

<!-- plan-path path=tests/testthat/test-matrix_adapter.R class=decision destination=@/tests/testthat/test-matrix_adapter.R order=23 reason=verify-matrix-adapter -->

Adopt adapter coverage only with order 22, including identifier alignment,
Cartesian-expansion prevention, output validation, and missing policies.

<!-- plan-path path=tests/testthat/test-step_result.R class=decision destination=@/tests/testthat/test-step_result.R order=24 reason=verify-structured-results -->

After orders 21 and 22, verify matrix/long structured results, assay metadata,
materialization, and `pb_eval()` behavior.
