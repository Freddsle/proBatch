# Project Donna workflows

## Goal of the document

This document describes the Donna workflows required to maintain proBatch as a specified, tested, documented, and releasable R/Bioconductor package.

## Scope

The scope of this specification includes permanent project-owned Donna workflows, their triggers, stages, verification responsibilities, composition, and completion outcomes.

Donna runtime internals, temporary session workflows, CI-provider configuration, and the detailed behavior of R package functions are out of scope.

## Dictionary

- `permanent workflow` — a version-controlled `./workflows/*.donna.md` artifact that implements a recurring project maintenance process.
- `planned workflow` — a required permanent workflow defined by this catalog whose workflow artifact has not yet been implemented.
- `verification gate` — a deterministic command or explicit review decision that must pass before a workflow can finish successfully.
- `repair loop` — a workflow transition from failed verification to focused agent repair and back to the same verification gate.
- `child workflow` — a workflow completed while another workflow remains active and resumes afterward.
- `reference-only migration entry` — an ordered source-commit record that preserves immutable history evidence and a reserved sequence slot but has no executable workflow artifact.

## General requirements

Permanent workflows MUST live under `./workflows/`, MUST end with `.donna.md`, and MUST be listed in this catalog.

Every permanent workflow MUST have a concise H1 summary, stable operation identifiers, an intentional start operation, at least one successful finish operation, and an explicit blocked or repair path for failures that require agent judgment.

Deterministic project commands MUST use `donna.lib.run_script`. Research, file edits, diagnosis, review, and decisions MUST use `donna.lib.request_action`.

Every automated failure that can be repaired within the workflow scope SHOULD enter a repair loop. A workflow MUST NOT report success merely because a command ran; its declared verification gates MUST pass.

Workflow scripts MUST run from the project root, use non-interactive commands, avoid network access unless the initiating request authorizes it, and preserve user changes outside the workflow scope.

A workflow MUST query Depmesh before changing governed artifacts when relevant relations are configured.

A workflow MUST NOT reset the Donna session. Starting a new session remains a separate developer-authorized action.

Every workflow artifact MUST pass `donna -p llm validate <artifact-id>` and SHOULD be rendered in view mode before its first execution and after control-flow changes.

Any workflow that drafts, proposes, creates, or amends a commit message MUST require compliance with `./specs/general/commits.md`. A workflow MUST NOT create or rewrite a commit unless the initiating request explicitly authorizes that Git action.

## Workflow reports

Permanent workflow output intended for later human action MUST be written as a Markdown report under `./workflows/reports/`. Reports MUST NOT use the `.donna.md` extension and MUST NOT be listed in the workflow catalog, because they are workflow output rather than workflows.

A finding that no remaining workflow resolves MUST be recorded as an entry in an open-questions report rather than left only in Donna session notes, which are cleared with the session.

`./workflows/reports/bec-ecoli-core-open-questions.md` owns the maintainer decisions left open by the BEC E. coli core migration. Each entry MUST identify the affected artifacts, the evidence, the effect, the workflows checked for coverage, and the candidate resolutions.

## Workflow catalog

### Bootstrap project specifications

- Artifact: `./workflows/bootstrap-project-specifications.donna.md`
- State: implemented
- Purpose: establish or revise the specification index, meta-specification, project agent guidance, file-relation model, Depmesh implementation, and workflow catalog.
- Triggers: initial specification-system setup; an explicit request to restructure the specification hierarchy; or a coordinated revision of the project relation and workflow model.
- Major stages: read reference specifications; inspect the R/Bioconductor package structure; author core specifications and agent guidance; specify and implement file relations; document required workflows; run deterministic deliverable checks.
- Expected verification: mandatory specification headings exist; `AGENTS.md` references the specification and package structure; Depmesh lists and resolves governance and test relations in both directions without accessing `./man/`; this workflow validates.
- Completion outcome: the specification system and dependency model are internally consistent, and planned workflows are clearly distinguished from implemented workflow artifacts.

### Migrate BEC E. coli core history

- Artifact: `./workflows/001_bec-ecoli-core-migration.donna.md`
- State: implemented
- Purpose: analyze the core-owned API and the stopped split reference, define the exact BEC history scope, generate numbered workflows for core-relevant commits, retain unrelated companion-only commits as ordered references, and coordinate source-first and split-reference migration with maintainer-controlled review and commit pauses.
- Triggers: the explicit migration of core-owned behavior from the local `BEC_ecoli_data` branch after baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32`.
- Major stages: verify pinned reference inputs; resolve the non-linear history boundary; inventory API ownership, commit effects, and split implementation evidence; classify each selected commit as executable or reference-only; create and catalog numbered child workflows; run executable children in order while skipping ordered references; pause twice per executable commit for maintainer analysis and commits; audit residual stopped-split changes; create structured residual reports; run deterministic governance checks.
- Expected verification: pinned commits and the API map match recorded identities; every selected commit appears exactly once in the manifest as either an executable child or a reference-only skip; every executable child appears exactly once in the workflow catalog, validates, and renders; reserved reference slots have no workflow artifacts and use only deterministic `reference-only` progress records; Depmesh returns exact workflow governance; focused source and test checks pass within each child; the maintainer explicitly resumes both review gates; final residual reports exist; generated-file-excluding `git diff --check` passes.
- Completion outcome: every selected commit has an explicit transferred, split-adjusted, reference-only, skipped, or blocked outcome; all package changes remain committed only by the maintainer; and remaining stopped-split candidates and API decisions are recorded with rationale and ordering.

At each maintainer-controlled review and commit pause in this migration family, the executing agent MUST inspect the exact in-scope unstaged diff and the contents of any in-scope untracked files without modifying the Git index. If that change set is nonempty, the agent MUST validate a complete copy-paste-ready commit message against `./specs/general/commits.md` and present it in chat in a fenced `text` block before waiting for the maintainer to resume. If one commit cannot accurately describe the change set, the agent MUST instead provide one validated copy-paste message per coherent commit and identify each intended file set. If the change set is empty, the agent MUST state that no commit message is needed and MUST NOT invent one. Presenting a message MUST NOT be treated as authorization to stage, create, or amend a commit, and unrelated pre-existing changes MUST remain excluded.

The BEC E. coli migration workflow family MUST use a zero-padded three-digit sequence followed by an underscore. `./workflows/001_bec-ecoli-core-migration.donna.md` owns the first number, and its generated foundation, per-commit, and residual workflows MUST receive subsequent unused numbers in execution order. Existing workflows outside this migration family MUST NOT be renamed solely to adopt this sequence.

Every permanent child generated by this migration MUST be added to this catalog with its exact artifact path before execution. Its state MUST remain `implementation in progress` until the artifact validates, renders correctly, and completes a representative execution.

For the selected post-synchronization history, a commit that changes only companion-owned implementation, tests, documentation, or adapters and has no core-owned editable artifact MUST remain in the ordered manifest as a reference-only record instead of receiving a permanent child workflow. A reference-only record MUST preserve its full source SHA, manifest ID, original sequence slot, classification, and concise change summary; its reserved prefix MUST NOT be reused or renumbered. It MUST NOT be cataloged as a workflow, executed, or given maintainer review pauses. When the parent reaches the entry, it MUST record only a deterministic `reference-only` progress outcome with no workflow and no source or split commit.

The reference-only commits and reserved prefixes are:

- `005_` — `65f70a46c4cf44e2717744aeafbf8acbe83b0378`: PRONE condition and batch metadata, backend compatibility, provider tests, and vignette.
- `018_` — `b8a262b4256966d60e1e8452ebde7a1bf471b4af`: omicsGMF batch-mean correction semantics and provider tests.
- `024_` — `4a06d99949114b2804b8d34492a288872fb611ed`: missing-policy forwarding through companion provider adapters and tests.
- `025_` — `6601232c69b44c507f2e9f63d256836e655f7973`: omicsGMF design-matrix column-name restoration.
- `026_` — `96d38eb7449e4d38c0d6a3fffd66e17f145f669f`: omicsGMF latent artifacts attached to provider outputs.
- `029_` — `0b3c3b55403f2bb7342a37236dd30f3d9b3544e9`: PRONE default input assay changed to `log2`.
- `030_` — `20e76c9b9b28c0ec98faa63c3f9382c1347301b9`: combined omicsGMF imputation/correction, aliases, and tests.

Commit `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57` at prefix `017_` MUST remain an executable child despite having no editable diff because its generated-only documentation change records core-related `handle_missing_values` evidence without allowing agents to inspect generated files.

#### BEC foundation core API

- Artifact: `./workflows/002_foundation_core_api.donna.md`
- State: implemented
- Purpose: migrate the 104 core-owned exports present at the synchronization merge while excluding companion providers and resolving effective correction load order.
- Triggers: the recommended post-synchronization scope selected by the parent migration.
- Major stages: inspect the pinned merge; port core behavior and tests; pause for source review; compare the stopped split; apply justified hardening; pause for split review; reverify.
- Expected verification: focused source/test checks, generated-file-excluding diff checks, and both post-review reverification gates pass.
- Completion outcome: the foundation core API has an explicit source and split result recorded by the parent.

#### BEC commit 42544a21f10c

- Artifact: `./workflows/003_42544a21f10c_pvca_impl_dedup.donna.md`
- State: implemented
- Purpose: review PVCA implementation consolidation and duplicate-definition removal.
- Triggers: ordered migration of source commit `42544a21f10ca6960d3e4c44d2833f764054d721`.
- Major stages: inspect the pinned PVCA hunks; port or record equivalence; review; compare split consolidation; review; reverify.
- Expected verification: PVCA sources/tests and symbol ownership pass focused checks without generated artifacts.
- Completion outcome: the PVCA commit is classified as transferred, adjusted, skipped, or blocked.

#### BEC commit db128668458a

- Artifact: `./workflows/004_db128668458a_prone_normalization_integration.donna.md`
- State: implemented
- Purpose: retain the core QFeatures link guard while excluding PRONE normalization providers.
- Triggers: ordered migration of source commit `db128668458a58ad31f66be5b6e39e2fedadbbe1`.
- Major stages: split core and companion hunks; verify the link behavior; review; compare registry/ownership evidence; review; reverify.
- Expected verification: core link tests pass and no PRONE implementation or dependency enters core.
- Completion outcome: the mixed commit has an explicit core and companion disposition.

#### BEC commit 4e4e9811503c

- Artifact: `./workflows/006_4e4e9811503c_grouped_na_heatmaps.donna.md`
- State: implemented
- Purpose: migrate grouped missingness heatmaps and their focused tests.
- Triggers: ordered migration of source commit `4e4e9811503cc7da1e35a8113a2fb8383fc5007b`.
- Major stages: inspect grouped heatmap APIs; port/test; review; compare split helpers; review; reverify.
- Expected verification: grouped aggregation, validation, and multi-assay tests pass.
- Completion outcome: grouped heatmap behavior has a reviewed core result.

#### BEC commit 12d4597e9e86

- Artifact: `./workflows/007_12d4597e9e86_violin_qualification.donna.md`
- State: implemented
- Purpose: qualify core split-violin calls while excluding companion variance-partition metadata.
- Triggers: ordered migration of source commit `12d4597e9e86ea2daddc77c4d236e692bbfa8fa9`.
- Major stages: split mixed hunks; port/test core calls; review; compare split source and dependencies; review; reverify.
- Expected verification: split-violin tests pass and `variancePartition` remains outside core.
- Completion outcome: the mixed metadata/source commit has an explicit disposition.

#### BEC commit 5747090a1de1

- Artifact: `./workflows/008_5747090a1de1_test_grouped_na_binarization.donna.md`
- State: implemented
- Purpose: migrate the grouped-heatmap binarization regression tests in source order.
- Triggers: ordered migration of source commit `5747090a1de1cdde780c6db7a84f988aedb9e8af`.
- Major stages: inspect test-only hunks; add or record equivalent coverage; review; compare split tests; review; reverify.
- Expected verification: threshold, ordering, and invalid-range cases pass after the following implementation is accounted for.
- Completion outcome: the test-only commit has a reviewed coverage result.

#### BEC commit 7feb6d5cea9e

- Artifact: `./workflows/009_7feb6d5cea9e_grouped_na_density.donna.md`
- State: implemented
- Purpose: migrate grouped density plots and grouped-heatmap binarization.
- Triggers: ordered migration of source commit `7feb6d5cea9ec92d378453ff65d672f231fa44b8`.
- Major stages: inspect source/tests; port; review; compare split public signatures; review; reverify.
- Expected verification: grouped density and binarization focused tests pass.
- Completion outcome: missingness visualization behavior has an explicit reviewed result.

#### BEC commit 0dbce8129444

- Artifact: `./workflows/010_0dbce8129444_density_color_scheme.donna.md`
- State: implemented
- Purpose: migrate grouped-density color schemes and editable documentation.
- Triggers: ordered migration of source commit `0dbce8129444c9079a560fb9bcd60b328b1f054c`.
- Major stages: inspect color hunks/tests; port; review; compare stopped-split helpers; review; reverify.
- Expected verification: named and palette color cases pass with the public grouped-density contract.
- Completion outcome: density color behavior has a reviewed result.

#### BEC commit 6d3c15e41905

- Artifact: `./workflows/011_6d3c15e41905_subset_and_bpca.donna.md`
- State: implemented
- Purpose: migrate sample subsetting and decide the optional Bayesian PCA surface.
- Triggers: ordered migration of source commit `6d3c15e419058d264b61ac43c27e95e0d635ca01`.
- Major stages: inspect container/PCA hunks; port; review; compare split behavior/dependencies; review; reverify.
- Expected verification: sample-subsetting and assay-preservation tests pass; the `pcaMethods` decision is explicit.
- Completion outcome: core subsetting and PCA behavior have reviewed results.

#### BEC commit 7db05faed18f

- Artifact: `./workflows/012_7db05faed18f_omicsgmf_all_na_guard.donna.md`
- State: implemented
- Purpose: retain core colData compatibility while excluding omicsGMF and PRONE provider behavior.
- Triggers: ordered migration of source commit `7db05faed18f75572d92553ae668288c326bda2a`.
- Major stages: separate mixed hunks; test core compatibility; review; compare ownership boundaries; review; reverify.
- Expected verification: factor/integer colData tests pass and provider implementations stay absent.
- Completion outcome: the mixed commit has an explicit core/companion result.

#### BEC commit fdafc6c8fe13

- Artifact: `./workflows/013_fdafc6c8fe13_centralize_all_na_handling.donna.md`
- State: implemented
- Purpose: assess shared all-NA helpers, correction wrappers, and normalization behavior without provider leakage.
- Triggers: ordered migration of source commit `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5`.
- Major stages: separate core/provider hunks; port/test; review; compare canonical split policies; review; reverify.
- Expected verification: retained correction, missing-data, normalization, and registry checks pass.
- Completion outcome: each mixed hunk has a reviewed destination or exclusion.

#### BEC commit 77a3365f0e13

- Artifact: `./workflows/014_77a3365f0e13_demote_all_na_helper_roxygen.donna.md`
- State: implemented
- Purpose: apply the internal-helper documentation change only if its prerequisite helpers are retained.
- Triggers: ordered migration of source commit `77a3365f0e13e9abca40a14e1898a351a1cdbc7f`.
- Major stages: inspect prerequisite state; change or skip source comments; review; compare split absence; review; reverify.
- Expected verification: editable R source parses and no generated documentation is accessed.
- Completion outcome: the documentation-only commit has a reviewed conditional result.

#### BEC commit 4b117a03f0c5

- Artifact: `./workflows/015_4b117a03f0c5_mask_groups_plot_na_intensity.donna.md`
- State: implemented
- Purpose: migrate group-failure masking and the NA-intensity plotting API.
- Triggers: ordered migration of source commit `4b117a03f0c5295f2253db26e13b945d7e3027b5`.
- Major stages: inspect source/tests; port; review; compare split implementations and coverage; review; reverify.
- Expected verification: grouped masking and NA-intensity focused tests pass.
- Completion outcome: both core behaviors have reviewed results.

#### BEC commit 4285c42f3167

- Artifact: `./workflows/016_4285c42f3167_clarify_removebatcheffect_missing_policy.donna.md`
- State: implemented
- Purpose: translate legacy removeBatchEffect missing-value behavior into the canonical core policy.
- Triggers: ordered migration of source commit `4285c42f31670d2f750dc8eb8c7ff1d0134a342d`.
- Major stages: inspect semantics; port or record supersession; review; compare split policy; review; reverify.
- Expected verification: correction, missing-policy, and matrix-adapter checks pass.
- Completion outcome: the core semantics are reviewed without mechanically restoring legacy booleans.

#### BEC commit 4540aca9182c

- Artifact: `./workflows/017_4540aca9182c_generated_missing_docs.donna.md`
- State: implemented
- Purpose: record a generated-only documentation commit without reading generated artifacts.
- Triggers: ordered migration of source commit `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57`.
- Major stages: verify the empty editable diff; review no-change; compare split provenance; review; reverify.
- Expected verification: no editable migration artifact exists and generated paths remain untouched.
- Completion outcome: the generated-only commit has a reviewed no-change result.

#### BEC commit d95dd736cb27

- Artifact: `./workflows/019_d95dd736cb27_harden_cv_dates_loess_mcombat.donna.md`
- State: implemented
- Purpose: migrate core CV, date, and LOESS safeguards while excluding mComBat.
- Triggers: ordered migration of source commit `d95dd736cb27d68fb2e21b20ab97d0f69826a663`.
- Major stages: split mixed hunks; port/test core safeguards; review; compare absent split behavior; review; reverify.
- Expected verification: CV, date-conversion, and nonlinear-fit regression checks pass.
- Completion outcome: core safeguards and companion exclusion have reviewed results.

#### BEC commit 6989542e4ece

- Artifact: `./workflows/020_6989542e4ece_harden_sample_alignment_consistency.donna.md`
- State: implemented
- Purpose: migrate independent sample-alignment, correlation, copy-safety, and consistency diagnostics.
- Triggers: ordered migration of source commit `6989542e4ece84f9fe210f3bf06b8254684e425e`.
- Major stages: separate mixed hunks; port/test; review; compare split equivalence and removals; review; reverify.
- Expected verification: correction, correlation, fit, and consistency checks pass for retained APIs.
- Completion outcome: each independent core hunk has a reviewed result.

#### BEC commit 428ad74b73fb

- Artifact: `./workflows/021_428ad74b73fb_cleanup_provider_control_flow_docs.donna.md`
- State: implemented
- Purpose: migrate core warning-capture and editable documentation while excluding provider and embedding cleanup.
- Triggers: ordered migration of source commit `428ad74b73fbc067de09d0a71c93395ae85eb51f`.
- Major stages: split mixed source/Roxygen hunks; port; review; compare split docs/control flow; review; reverify.
- Expected verification: retained R sources and focused tests pass without generated documentation.
- Completion outcome: core documentation/control-flow changes and provider exclusions are explicit.

#### BEC commit 72d11d1f7f92

- Artifact: `./workflows/022_72d11d1f7f92_version_variancepartition.donna.md`
- State: implemented
- Purpose: review release version metadata and exclude the companion variancePartition dependency.
- Triggers: ordered migration of source commit `72d11d1f7f92a10cb9a5244e68f9c19e830197b9`.
- Major stages: inspect metadata; record current-policy result; review; compare split metadata; review; reverify.
- Expected verification: `DESCRIPTION` remains valid and contains no companion-only dependency.
- Completion outcome: both metadata hunks have explicit reviewed decisions.

#### BEC commit 8088de1701e0

- Artifact: `./workflows/023_8088de1701e0_ignore_docker_readme.donna.md`
- State: implemented
- Purpose: review repository-local Docker ignore metadata.
- Triggers: ordered migration of source commit `8088de1701e0908cca25b978105f8d6c7bfccc20`.
- Major stages: inspect `.gitignore`; change or skip; review; compare split metadata; review; reverify.
- Expected verification: ignore syntax and generated-file-excluding diff checks pass.
- Completion outcome: the metadata-only commit has a reviewed result.

#### BEC commit b60e2b169afe

- Artifact: `./workflows/027_b60e2b169afe_store_explicit_final_log_assay.donna.md`
- State: implemented
- Purpose: migrate or record equivalent explicit final log-assay storage.
- Triggers: ordered migration of source commit `b60e2b169afe45934c18eee1c83469e7e5fed33f`.
- Major stages: inspect registry/container hunks; port/test; review; compare strengthened split behavior; review; reverify.
- Expected verification: explicit final log assay, naming, and lineage checks pass.
- Completion outcome: the core behavior has a reviewed transferred or equivalent result.

#### BEC commit 5edcc4b89f5c

- Artifact: `./workflows/028_5edcc4b89f5c_report_inplace_link_effects.donna.md`
- State: implemented
- Purpose: preserve correct assay-link behavior without replaying obsolete warnings.
- Triggers: ordered migration of source commit `5edcc4b89f5cbfba07e1ca1057ab475a56501202`.
- Major stages: inspect link-warning hunks; port or skip; review; compare evolved split behavior; review; reverify.
- Expected verification: missing-filter and assay-link tests pass.
- Completion outcome: each filter's current link behavior and warning decision are reviewed.

#### BEC commit e2bb18547c73

- Artifact: `./workflows/031_e2bb18547c73_drop_unused_model_levels.donna.md`
- State: implemented
- Purpose: migrate unused-factor-level handling for core correction while excluding BERT modeling.
- Triggers: ordered migration of source commit `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.
- Major stages: separate mixed core/BERT hunks; port/test; review; compare split gap; review; reverify.
- Expected verification: ComBat/removeBatchEffect unused-level regressions pass and BERT stays outside core.
- Completion outcome: the tip commit has explicit core and companion dispositions.

#### Residual stopped-split review

- Artifact: `./workflows/032_residual_split_review.donna.md`
- State: implemented
- Purpose: classify the immutable 99-path residual split delta and create a structured remaining-change plan.
- Triggers: completion of the foundation and 22 executable post-synchronization commit children, with all seven reference-only commits recorded as deterministic skips.
- Major stages: inspect pinned split evidence; draft two reports; verify exact path/API markers; repair; pause for report review; reverify.
- Expected verification: 99 ordered path markers reproduce the pinned checksum, 26 companion API markers are exact, and every actionable path has one plan marker.
- Completion outcome: residual recommendations, exclusions, and API decisions are reviewed and recorded for the parent.

### Verify specifications

- Artifact: `./workflows/verify-specifications.donna.md`
- State: planned
- Purpose: verify structural and semantic consistency across all project specifications without changing implementation behavior.
- Triggers: any change under `./specs/`; a specification path change; or a pre-release documentation review.
- Major stages: compare the filesystem with `./specs/intro.md`; check mandatory headings and path style; review RFC 2119 statements for clarity and conflicts; query governance relations; repair discrepancies; repeat checks.
- Expected verification: every specification is indexed exactly once; every specification contains `Goal of the document` and `Scope`; all except the meta-specification are governed by it; `git diff --check` passes; Depmesh queries match the current paths.
- Completion outcome: specifications are structurally valid, indexed, non-conflicting, and connected to governed artifacts.

### Verify file relations

- Artifact: `./workflows/verify-file-relations.donna.md`
- State: planned
- Purpose: validate the complete operational contract in `./specs/behavior/files_relations.md`.
- Triggers: changes to `./depmesh.toml`, `./bin/depemesh/`, specification paths, R/test filename conventions, documentation sources under `./R/`, or permanent workflow paths.
- Major stages: list relation definitions; query conventional and exceptional test mappings in both directions; query governance in both directions; verify that no rule or helper accesses `./man/`; repair rules or helpers; rerun all checks.
- Expected verification: every required relation is listed; representative and exceptional mappings return exact existing artifacts; no relation emits paths outside the project or reads, enumerates, returns, or checks files under `./man/`.
- Completion outcome: Depmesh provides deterministic, bidirectional impact discovery consistent with the file-relation specification.

### Run focused package tests

- Artifact: `./workflows/run-tests.donna.md`
- State: planned
- Purpose: run the package’s testthat suite and guide focused repair without the broader cost of a full package check.
- Triggers: changes to `./R/`, `./tests/`, test fixtures, or package data used by tests; or an explicit request for fast regression verification.
- Major stages: query source/test relations; select focused tests when the change scope is known; run focused tests; run the complete testthat suite; diagnose failures; repair in scope; repeat until green.
- Expected verification: affected focused tests pass and the complete testthat suite exits successfully in a clean R process.
- Completion outcome: testthat regression coverage passes, with skipped tests and environment limitations reported explicitly.

### Review documentation sources

- Artifact: `./workflows/review-documentation-sources.donna.md`
- State: planned
- Purpose: review and improve public API documentation at its editable source under `./R/` without reading or generating `./man/`.
- Triggers: changes to Roxygen2 comments, exports, S3 or S4 registration, public signatures, documented datasets, or documentation-source conventions.
- Major stages: inspect affected R files and public interfaces; edit Roxygen2 comments under `./R/`; review source-level completeness and consistency; confirm no file under `./man/` was accessed or changed; hand generation off to the maintainer.
- Expected verification: editable Roxygen2 sources are internally consistent with the affected R interfaces, and the workflow has not read, generated, validated, compared, or modified `./man/`.
- Completion outcome: documentation sources under `./R/` are ready for the maintainer to regenerate manually with devtools; generated output remains uninspected.

### Build package vignettes

- Artifact: `./workflows/build-vignettes.donna.md`
- State: planned
- Purpose: verify executable long-form documentation and its supporting assets.
- Triggers: changes under `./vignettes/`, changes to APIs or example data used by vignettes, documentation dependency changes, or release preparation.
- Major stages: inspect vignette dependencies and affected package APIs; build vignettes in a clean process; capture rendering failures; repair code, prose, or assets; rebuild; confirm generated caches and HTML are not accidentally tracked.
- Expected verification: every declared vignette builds without error using package-declared dependencies, references and assets resolve, and only intended source artifacts remain in the working tree.
- Completion outcome: package vignettes are executable, current with the public API, and ready for package building.

### Check R and Bioconductor package

- Artifact: `./workflows/check-package.donna.md`
- State: planned
- Purpose: perform the authoritative local quality gate for an R/Bioconductor package.
- Triggers: completion of a nontrivial implementation change; changes to metadata, exports, compiled or generated artifacts, dependencies, examples, tests, or vignettes; or release preparation.
- Major stages: inspect package metadata and dependency impact; review documentation sources when required; run focused agent-safe source and test checks; ask the maintainer to regenerate documentation and run any package-wide R or Bioconductor checks that inspect `./man/`; repair only non-`man/` issues in scope; repeat the applicable checks.
- Expected verification: agent-safe checks pass; the maintainer confirms completion of package-wide checks; any `man/` diagnostics are reported and left entirely to the maintainer; the source tree contains no unintended check/build artifacts.
- Completion outcome: agent-owned checks have passed and maintainer-owned generation and package-wide verification have been explicitly handed off or confirmed without agent inspection of `./man/`.

### Prepare package release

- Artifact: `./workflows/prepare-release.donna.md`
- State: planned
- Purpose: coordinate all artifacts and checks needed for a deliberate proBatch release without publishing it.
- Triggers: an explicit release-preparation request with a target version and release scope.
- Major stages: verify version and dependency metadata in `./DESCRIPTION`; update `./NEWS` and editable user-facing documentation; run specification and file-relation verification; review documentation sources; run tests and vignette builds; obtain maintainer confirmation for documentation generation and package-wide checks; review the final editable-source diff.
- Expected verification: every agent-owned child workflow completes; version records agree; release notes cover user-visible changes; maintainer-owned generation and checks are confirmed; no runtime, cache, or check artifacts are included.
- Completion outcome: the repository and verified source archive are ready for a separately authorized release or submission action; the workflow MUST NOT push, tag, upload, or submit by itself.

## Workflow composition

`./workflows/check-package.donna.md` SHOULD invoke the documentation-source, test, and vignette workflows as child workflows when their checks cannot be performed equivalently by its deterministic agent-safe stages.

`./workflows/prepare-release.donna.md` MUST coordinate specification verification, file-relation verification, documentation-source review, tests, vignettes, and maintainer-owned generation and package checks. It SHOULD invoke the specialized permanent workflows rather than duplicating their repair logic.

The specialized workflows MUST remain independently runnable so maintainers can obtain focused feedback without executing the complete release pipeline.

## Implementation order

Planned workflows SHOULD be implemented in the following order:

1. `./workflows/verify-file-relations.donna.md`
2. `./workflows/verify-specifications.donna.md`
3. `./workflows/run-tests.donna.md`
4. `./workflows/review-documentation-sources.donna.md`
5. `./workflows/build-vignettes.donna.md`
6. `./workflows/check-package.donna.md`
7. `./workflows/prepare-release.donna.md`

A planned workflow MUST NOT be presented as executable until its artifact exists, validates, renders correctly, and has completed at least one representative execution.

When a workflow is implemented, renamed, superseded, or removed, this catalog, affected Depmesh rules, and `./AGENTS.md` MUST be updated in the same change.
