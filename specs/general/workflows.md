# Project Donna workflows

## Goal of the document

This document describes the Donna workflows required to maintain proBatch as a specified, tested, documented, standalone, and releasable R/Bioconductor package.

## Scope

The scope of this specification includes permanent project-owned Donna workflows, their triggers, stages, verification responsibilities, composition, completion outcomes, and preservation of completed workflow history.

Donna runtime internals, temporary session workflows, CI-provider configuration, and the detailed behavior of R package functions are out of scope.

## Dictionary

- `permanent workflow` — a version-controlled `./workflows/*.donna.md` artifact that implements a recurring or bounded project process.
- `planned workflow` — a required permanent workflow defined by this catalog whose workflow artifact has not yet been implemented.
- `verification gate` — a deterministic command or explicit review decision that must pass before a workflow can finish successfully.
- `repair loop` — a workflow transition from failed verification to focused agent repair and back to the same verification gate.
- `child workflow` — a workflow completed while another workflow remains active and resumes afterward.
- `historical workflow record` — a non-executable Markdown report that consolidates the immutable inputs, operations, outcomes, verification, and superseding decisions of a removed completed workflow family.

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

## Workflow reports and history

Permanent workflow output intended for later human action MUST be written as a Markdown report under `./workflows/reports/`. Reports MUST NOT use the `.donna.md` extension and MUST NOT be listed in the workflow catalog, because they are workflow output rather than workflows.

A finding that no remaining workflow resolves MUST be recorded as an entry in an open-questions report rather than left only in Donna session notes, which are cleared with the session.

A completed workflow family MAY be removed after one historical workflow record preserves its immutable inputs, ordered operations, outcomes, destination changes, verification evidence, exclusions, and later superseding decisions. The historical record MUST state that it is non-executable and MUST NOT be used to reconstruct removed workflow state.

`./workflows/reports/bec-ecoli-core-migration-history.md` is the sole historical workflow record for the completed BEC E. coli core migration. It preserves the migration evidence formerly distributed across the numbered workflow family, manifest, progress ledger, and residual reports.

## Workflow catalog

### Bootstrap project specifications

- Artifact: `./workflows/bootstrap-project-specifications.donna.md`
- State: implemented
- Purpose: establish or revise the specification index, meta-specification, project agent guidance, file-relation model, Depmesh implementation, and workflow catalog.
- Triggers: initial specification-system setup; an explicit request to restructure the specification hierarchy; or a coordinated revision of the project relation and workflow model.
- Major stages: read reference specifications; inspect the R/Bioconductor package structure; author core specifications and agent guidance; specify and implement file relations; document required workflows; run deterministic deliverable checks.
- Expected verification: mandatory specification headings exist; `./AGENTS.md` references the specification and package structure; Depmesh lists and resolves governance and test relations in both directions without accessing `./man/`; this workflow validates.
- Completion outcome: the specification system and dependency model are internally consistent, and planned workflows are clearly distinguished from implemented workflow artifacts.

### Establish standalone core baseline

- Artifact: `./workflows/establish-standalone-core-baseline.donna.md`
- State: implementation in progress
- Purpose: turn the migrated package into a coherent standalone `proBatch` baseline that is safe for `proBatchBench` to import, while restoring t-SNE and UMAP diagnostics and avoiding a whole-file replay of the rejected split rewrite.
- Triggers: the maintainer's explicit request to close the core migration findings, restore `plot_TSNE()` and `plot_UMAP()` from original proBatch, and establish the one-way Core-to-Bench ownership boundary.
- Major stages: verify pinned source and downstream evidence; specify the standalone ownership and extension contracts; repair package hygiene and core invariants; implement only the downstream extension surface justified by actual Bench use; restore t-SNE and UMAP with focused tests; synchronize dependencies, relations, Roxygen sources, vignettes, and release notes; run agent-safe verification; obtain maintainer confirmation for generated documentation and package-wide checks; verify the installed import surface.
- Expected verification: the behavior specification and Depmesh governance are consistent; focused and complete testthat runs pass; `./DESCRIPTION`, Pixi metadata, dependency usage, and build exclusions agree; t-SNE and UMAP work for matrix and `ProBatchFeatures` inputs with guarded optional rendering; Core has no dependency on Bench or provider engines; Donna and Depmesh validation pass; the maintainer confirms generated documentation, vignette, package, Bioconductor, ComBat, and quantile-normalization checks; the loaded namespace exposes the agreed import surface without inspecting `./man/`.
- Completion outcome: `proBatch` is a verified standalone baseline with a minimal documented extension contract, Core-owned PCA/t-SNE/UMAP diagnostics, and an explicit handoff for removing duplicated implementations from `proBatchBench`.

### Verify specifications

- Artifact: `./workflows/verify-specifications.donna.md`
- State: implemented
- Purpose: verify structural and semantic consistency across all project specifications without changing implementation behavior.
- Triggers: any change under `./specs/`; a specification path change; or a pre-release documentation review.
- Major stages: capture the initial repository scope excluding `./man/`; compare the specification filesystem with `./specs/intro.md`; check mandatory headings and path style; review RFC 2119 statements for clarity and conflicts; query governance relations; repair discrepancies; snapshot the reviewed authorized edit set; repeat protected-scope checks immediately before completion.
- Expected verification: every specification is indexed exactly once; every specification contains `Goal of the document` and `Scope`; all except the meta-specification are governed by it; `git diff --check` passes; Depmesh queries match the current paths; HEAD, staged entries, assume-unchanged, skip-worktree, fsmonitor-valid, intent-to-add, and resolve-undo state outside `./man/` remain unchanged; tracked and non-ignored untracked artifacts outside the authorized edit set remain at their initial state; and the authorized edit set remains at its reviewed state through completion.
- Completion outcome: specifications are structurally valid, indexed, non-conflicting, and connected to governed artifacts.

### Verify file relations

- Artifact: `./workflows/verify-file-relations.donna.md`
- State: implemented
- Purpose: validate the complete operational contract in `./specs/behavior/files_relations.md`.
- Triggers: changes to `./depmesh.toml`, `./bin/depemesh/`, specification paths, R/test filename conventions, documentation sources under `./R/`, or permanent workflow paths.
- Major stages: list relation definitions; query conventional and exceptional test mappings in both directions; query governance in both directions; verify that no rule or helper accesses `./man/`; repair rules or helpers; rerun all checks.
- Expected verification: every required relation is listed; representative and exceptional mappings return exact existing artifacts; no relation emits paths outside the project or reads, enumerates, returns, or checks files under `./man/`.
- Completion outcome: Depmesh provides deterministic, bidirectional impact discovery consistent with the file-relation specification.

### Run focused package tests

- Artifact: `./workflows/run-tests.donna.md`
- State: implemented
- Purpose: run the package's testthat suite and guide focused repair without the broader cost of a full package check.
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
- Major stages: inspect package metadata and dependency impact; review documentation sources when required; run focused agent-safe source and test checks; ask the maintainer to regenerate documentation and run any package-wide R or Bioconductor checks that inspect `./man/`; repair only in-scope issues unrelated to `./man/`; repeat the applicable checks.
- Expected verification: agent-safe checks pass; the maintainer confirms completion of package-wide checks; any `./man/` diagnostics are reported and left entirely to the maintainer; the source tree contains no unintended check/build artifacts.
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
