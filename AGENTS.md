# Agent Instructions

## Specifications and project structure

The authoritative specification index is `specs/intro.md`. Before changing a
project artifact, read the index, the specifications it identifies as relevant,
and the dependencies reported by Depmesh. Treat normative statements in the
specifications as project requirements. If instructions and specifications
conflict, stop and report the conflict instead of silently choosing one.

All specification documents must follow `specs/meta/general.md`. Add, move,
rename, or remove a specification only when `specs/intro.md` and affected
Depmesh rules are updated in the same change.

This repository is an R-language Bioconductor package. Preserve these ownership
boundaries:

- `DESCRIPTION` defines package metadata, R/Bioconductor dependencies, and the
  supported R version.
- `R/` contains package implementation and the Roxygen2 comments that are the
  source of package documentation.
- `man/` is opaque, maintainer-owned generated output. Never read, edit,
  regenerate, lint, validate, diff-review, or use files under `man/` for
  dependency discovery. Make documentation changes only in Roxygen2 comments
  under `R/`; the maintainer manually regenerates `man/` with devtools.
- `NAMESPACE` is generated with the documentation. Change its Roxygen2 source
  under `R/` and leave regeneration to the maintainer because documentation
  generation also writes `man/`.
- `tests/testthat/` contains focused testthat tests; `tests/testthat.R` is the
  package test entry point.
- `vignettes/` contains executable package guides and their supporting assets.
- `data/` contains packaged example data, and `inst/` contains installed package
  resources.
- `specs/` contains project requirements, while `workflows/` contains
  project-owned Donna workflows.
- `bin/depemesh/` contains project-local helpers used only for dependency
  discovery that cannot be expressed directly in `depmesh.toml`.

For R changes, preserve public API documentation and registration conventions
in their `R/` sources, add or update corresponding testthat coverage, and verify
at the narrowest useful level before running broader package checks. Do not run
`devtools::document()`, `roxygen2::roxygenise()`, or any command that reads,
rewrites, or validates `man/`, including as part of a broader package check.
Leave any such package-wide command to the maintainer. If the maintainer reports
a `man/` problem, defer it without investigating or correcting the generated
file.

## Commit messages

Every commit message an agent drafts, proposes, creates, or amends must follow
`specs/general/commits.md` and Conventional Commits 1.0.0:

```text
<type>[optional scope][optional !]: <description>
```

Use only the types allowed by the specification. Before generating a message,
inspect the exact change set it will describe: at a maintainer-controlled
workflow pause, inspect the exact in-scope unstaged diff and any in-scope
untracked files without modifying the index; otherwise inspect the exact staged
diff or the explicit change summary. Choose the type from the change’s effect,
and validate the header, scope, length, description, and breaking-change
markers. If one message cannot accurately describe the inspected changes,
report that the commit should be split. Never invent issue references or
attribution, and never create or rewrite a commit without user authorization.

## Donna workflows

This project uses Donna as its workflow engine. The project configuration is
`donna.toml`; project-owned workflows belong in `workflows/`, and temporary
workflow/session artifacts belong in `.session/donna/`.

`specs/general/workflows.md` is the authoritative workflow catalog. Read it
before selecting or creating a project workflow. Catalog entries marked
`planned` are requirements, not executable workflows; do not claim or attempt
to run them until their artifacts exist and validate.

Before using Donna in a session:

1. Read the built-in instructions with `donna skill` (or
   `donna -p llm skill usage`).
2. Use the `llm` output protocol for agent-driven commands.
3. Run `donna -p llm status` before deciding whether to start or continue work.

Follow these rules:

- Address pending Donna action requests before unrelated work.
- If new work is requested while Donna has pending work units or action
  requests, ask whether to continue the current session or start a new one.
- Run workflows only when explicitly requested by the user, these project
  instructions, or Donna itself.
- Keep every permanent workflow under `workflows/` synchronized with the
  workflow catalog and its Depmesh governance relations.
- Treat completed-workflow history reports as non-executable evidence. Do not
  reconstruct or run removed workflow artifacts from a historical report.
- Keep external source and downstream repositories read-only unless the user
  explicitly places changes to them in scope.
- Copy action-request IDs and next-operation IDs exactly from Donna's output;
  never invent them.
- Do not run `donna -p llm new-session` unless resetting the active session is
  understood and explicitly allowed.
- Use project-root-anchored workflow IDs such as
  `@/workflows/example.donna.md`.

Common commands:

```sh
donna -p llm status
donna -p llm list
donna -p llm run @/workflows/example.donna.md
donna -p llm continue
donna -p llm validate --all
```

## Dependency discovery

This project uses Depmesh to discover specification governance and R
source/test coverage. The normative relation model is
`specs/behavior/files_relations.md`, and its implementation is `depmesh.toml`.

Before changing an artifact:

1. Read the built-in instructions with `depmesh skill` (or
   `depmesh -p llm skill usage`).
2. Run `depmesh -p llm relations` before relying on relation names.
3. Query all relations for the artifact with `depmesh -p llm dependencies` and
   inspect the returned files before editing. Use relation filters only when the
   task boundary already justifies a narrower query.

Use the `llm` output protocol for agent-driven commands and prefer
project-root-anchored artifact IDs such as `@/R/normalize.R`. Combine relation
filters when only a subset is relevant.

Common commands:

```sh
depmesh -p llm relations
depmesh -p llm dependencies @/R/normalize.R
depmesh -p llm dependencies --relation governed_by @/depmesh.toml
depmesh -p llm dependencies --relation tested_by @/R/normalize.R
depmesh -p llm dependencies --relation tests @/tests/testthat/test-normalize.R
```
