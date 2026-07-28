# Project file relations

## Goal of the document

This document describes how dependency relations between proBatch project artifacts are named, represented, discovered, queried, maintained, and verified.

## Scope

The scope of this specification includes relations between specifications, agent instructions, R source files, testthat tests, Donna workflows, Donna and Depmesh configuration, and dependency-discovery helpers.

Generated files under `./man/`, R function call graphs, runtime package dependencies, data provenance inside packaged datasets, and informal textual references that cannot be discovered reliably are out of scope.

## Dictionary

- `artifact` — a file in the project that can be queried through Depmesh with a root-anchored identifier.
- `artifact identifier` — a project-root-anchored path in the form `@/path/to/file`.
- `relation` — a named, directional connection from a queried artifact to one or more dependency artifacts.
- `forward relation` — the relation whose direction follows the primary maintenance question, such as an R source file being `tested_by` a test.
- `reverse relation` — the separately named inverse lookup, such as a test that `tests` an R source file.
- `governing specification` — a specification containing normative requirements that apply to an artifact.

## General principles

The normative relation model MUST be implemented in `./depmesh.toml`.

Every project file shown in Depmesh output MUST use a root-anchored artifact identifier beginning with `@/`.

Relation identifiers MUST use lowercase snake case. Their descriptions MUST state what is returned for the queried artifact rather than restating the relation name.

Relations MUST be directional. A reverse relation MUST be declared when maintainers need to navigate from either side of the relationship.

Rules MUST prefer scalable path patterns over enumerated files when the relationship follows a stable project convention. Exceptions to a convention MUST use explicit rules and MUST NOT weaken the conventional rule.

Rules MUST return only project-local artifacts. A relation command MUST NOT emit installed package paths, temporary files, URLs, or artifacts outside the repository.

An empty relation result means that no matching relation is currently known. It MUST NOT be interpreted as proof that an artifact has no conceptual dependencies.

## Required relation vocabulary

### `governed_by` and `governs`

`governed_by` MUST return specifications whose normative requirements apply to the queried artifact.

`governs` MUST return artifacts to which the queried specification applies.

Every specification except `./specs/meta/general.md` MUST be governed by `./specs/meta/general.md`.

`./AGENTS.md` MUST be governed by all current specifications because it directs agents to honor the specification system.

`./depmesh.toml`, dependency-discovery helpers, and artifact families whose relations are defined here MUST be governed by `./specs/behavior/files_relations.md`.

`./donna.toml`, project-owned `./workflows/*.donna.md` files, and workflow report artifacts under `./workflows/reports/` MUST be governed by `./specs/general/workflows.md`.

Reverse governance rules MUST expose the same relationships through `governs`.

### `tested_by` and `tests`

`tested_by` MUST return testthat files that directly or integratively verify the queried R source file.

`tests` MUST return R source files directly or integratively verified by the queried testthat file.

The conventional mapping MUST pair `./R/<module>.R` with `./tests/testthat/test-<module>.R` when both files exist.

Filename exceptions MUST be declared explicitly in both directions. A many-to-one integration test MUST list every source file whose coordinated behavior it intentionally verifies.

The presence of a test relation MUST NOT imply complete behavioral coverage. It records ownership and expected impact, not a coverage percentage.

## Artifact-family organization

### Specifications

Specifications live under `./specs/` and MUST follow `./specs/meta/general.md`.

The specification index at `./specs/intro.md` MUST enumerate every specification. Adding, moving, renaming, or removing a specification MUST update the index and affected governance rules in the same change.

### R sources and tests

Package implementation belongs under `./R/`. Focused tests belong under `./tests/testthat/` and SHOULD follow the conventional source/test stem mapping.

When a test intentionally spans multiple source files, its reverse `tests` result MUST enumerate the stable source ownership set. Incidental helper calls made during test execution MUST NOT be added as relations.

When a source or test filename cannot follow the convention, the exception MUST be recorded next to the conventional rules in `./depmesh.toml`.

### Maintainer-generated manuals

Public API documentation MUST originate in Roxygen2 comments under `./R/`.

Files under `./man/` MUST be excluded from all Depmesh inputs, outputs, helpers, and validation queries. Depmesh MUST NOT read, enumerate, validate, or infer relations from them.

The maintainer alone generates `./man/` with devtools. Dependency discovery MUST stop at the editable R and Roxygen2 sources under `./R/`.

### Donna workflows

Permanent project workflows MUST live under `./workflows/` and end with `.donna.md`. Temporary workflows MUST live under the Donna session directory configured in `./donna.toml`.

Only permanent workflows MUST be returned by project governance relations. Session artifacts MUST remain ignored runtime state and MUST NOT appear as project dependencies.

### Workflow reports

Permanent Markdown artifacts that a workflow produces for later human action MUST live under `./workflows/reports/` and MUST NOT use the `.donna.md` extension, so that Donna does not discover them as workflows.

A workflow report MUST record findings, decisions, or open questions rather than executable process steps. Findings that no remaining workflow resolves MUST be recorded in a workflow report instead of being left only in session notes.

Workflow reports MUST be governed by `./specs/general/workflows.md` and MUST NOT require a workflow catalog entry, because they are workflow output rather than workflows.

### Dependency-discovery helpers

Project-local Depmesh helpers MUST live under `./bin/depemesh/`.

Helpers MUST be deterministic, read-only, and executable from the project root. They MUST NOT access the network, modify project files, install packages, or depend on user-specific state.

A helper used as a Depmesh command source MUST write one artifact identifier per standard-output line, write diagnostics to standard error, emit no unrelated text, and exit nonzero for invalid invocation or malformed required input.

Helpers SHOULD use base R when they inspect R package metadata and no external package is required.

Project-only specifications, workflows, agent instructions, dependency configuration, session state, and helper scripts MUST be excluded from the built R package through `./.Rbuildignore`.

## Depmesh configuration structure

`./depmesh.toml` MUST declare relation definitions before rules.

Paired forward and reverse relations SHOULD be adjacent in the relation list. Rules SHOULD be grouped by relation family in this order:

1. specification governance;
2. source/test ownership;
3. additional stable relation families introduced by later specifications.

Within a relation family, the scalable conventional rule SHOULD appear before explicit exceptions.

Command rules MUST quote artifact arguments and MUST call helpers with paths relative to the project root.

## Query requirements

Agents MUST run `depmesh -p llm relations` before relying on relation identifiers in a work session.

Before editing an artifact, an agent MUST query all configured relations for that artifact unless a narrower query is justified by an already understood task boundary.

Agents SHOULD query reverse relations when changing a specification, test, generated document, or workflow catalog because those artifacts commonly affect multiple dependents.

Multiple input artifacts MAY be queried together for an impact-set overview. Separate queries MUST be used when attribution to each input artifact matters.

## Validation

Every change to `./depmesh.toml` or `./bin/depemesh/` MUST validate the relation list and representative queries in both directions.

Validation MUST include:

- a conventionally named R source and test pair;
- every explicit filename exception affected by the change;
- a specification-to-governed-artifacts reverse query;
- a governed project configuration or workflow artifact.

Validation MUST confirm that no relation rule or helper reads, enumerates, returns, or otherwise checks files under `./man/`.

Depmesh relation checks MUST be included in Donna workflows that modify specifications, package structure, documentation sources under `./R/`, or dependency configuration.
