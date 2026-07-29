# General specification requirements

## Goal of the document

This document describes the common structure, language, abstraction, and maintenance conventions for specifications in proBatch.

## Scope

The scope of this specification is limited to requirements for Markdown documents under `./specs/`.

Package behavior, implementation details, user documentation content, and workflow-specific requirements are out of scope except where they affect how specifications are authored or maintained.

## Dictionary

- `specification` — a Markdown document under `./specs/` that defines requirements, behavior, architecture, terminology, artifact relations, or maintenance rules for the project.
- `specification index` — `./specs/intro.md`, the authoritative list of specification directories and documents.
- `top-level section` — a section introduced by an `h2` Markdown heading.
- `nested section` — a section introduced by an `h3` or deeper Markdown heading.
- `normative statement` — a requirement expressed with the RFC 2119 key words `MUST`, `MUST NOT`, `SHOULD`, `SHOULD NOT`, or `MAY`.

## Document structure

A specification MUST contain exactly one `h1` heading that names the specification, and that name SHOULD be unique across the project.

Top-level information SHOULD be organized with `h2` headings. Nested sections MAY use `h3` and deeper headings when they make a larger topic easier to navigate.

Every specification MUST begin with the following top-level sections in this order:

1. `Goal of the document`
2. `Scope`
3. `Dictionary`, when specification-specific terminology is required

The `Goal of the document` section MUST describe the content and purpose of the specification. It MUST NOT define requirements for the document itself by saying that the document MUST list, define, or describe something.

The `Scope` section MUST describe the boundaries of the specification. It SHOULD identify important out-of-scope topics when doing so prevents ambiguity, and it MUST NOT redirect requirements to other specifications.

The `Dictionary` section SHOULD define only terms that are specific to the current specification. A term shared by multiple specifications SHOULD be defined in a dedicated project dictionary when one exists.

## Normative language and style

Specifications MUST use Markdown syntax.

Normative statements MUST use the key words defined by [RFC 2119](https://datatracker.ietf.org/doc/html/rfc2119).

Requirements MUST be testable or reviewable. A specification MUST avoid vague qualifiers such as “usually,” “where practical,” or “as appropriate” unless it also defines the condition that determines compliance.

Specifications MUST NOT insert hard line breaks merely to fit a fixed column width; paragraphs MUST use as many characters per source line as needed to express the idea clearly.

Long enumerations SHOULD use Markdown lists. Examples SHOULD be short, representative, and explicitly distinguishable from normative requirements.

Project-root-relative paths MUST use the `./path/to/artifact` form in prose. Depmesh artifact identifiers MUST use the root-anchored `@/path/to/artifact` form when they are shown as query inputs or outputs.

## Abstraction level

Specifications MUST describe behavior, architecture, constraints, terminology, and compatibility contracts at the highest level that remains precise enough to guide implementation and review.

Specifications SHOULD define:

- externally visible behavior and data contracts;
- stable ownership boundaries between package artifact families;
- invariants that must hold across implementations;
- intentional R or Bioconductor technology choices;
- examples that clarify a requirement without enumerating the entire current implementation.

Specifications MUST NOT define incidental implementation details such as private helper names, temporary algorithms, local variable names, or exact code locations unless those names or locations are stable project contracts.

Specifications MAY name concrete files, directories, R symbols, commands, configuration fields, or formats when those names form a stable interface or ownership boundary.

When a requirement can be expressed as either an implementation detail or a general architectural rule, the specification MUST use the general rule.

## R and Bioconductor context

Specifications that govern package code MUST account for the package metadata in `./DESCRIPTION`, the R source and Roxygen2 documentation comments under `./R/`, testthat tests under `./tests/`, executable vignettes under `./vignettes/`, and generated `./NAMESPACE`.

Files under `./man/` are opaque, maintainer-owned generated output. The maintainer manually generates them with devtools from R files and Roxygen2 comments under `./R/`.

Agents and agent-driven workflows MUST NOT read, edit, regenerate, lint, validate, compare, review, or use files under `./man/` for dependency discovery. They MUST NOT run `devtools::document()`, `roxygen2::roxygenise()`, or another command that may read or write `./man/`.

Documentation changes MUST be made only in R files and Roxygen2 comments under `./R/`. Generation of `./man/` and the associated `./NAMESPACE` output MUST be left to the maintainer.

An agent MUST NOT run a broader package command that reads or validates `./man/`; that command MUST be left to the maintainer. If a maintainer-run command reports a problem involving `./man/`, an agent MUST defer it without inspecting or correcting the generated file.

Requirements for public R APIs SHOULD describe observable function, S3, S4, or Bioconductor container behavior rather than internal implementation mechanics.

Requirements that affect package dependencies, exported interfaces, documentation sources, examples, tests, or vignettes MUST identify the affected editable artifact families so dependency discovery and verification can cover them without accessing `./man/`.

## Specification relationships

Every specification MUST be listed in `./specs/intro.md`. When a specification file or directory is added, moved, renamed, or removed, the index MUST be updated in the same change.

A specification SHOULD own one cohesive concern. Requirements MUST NOT be duplicated across specifications; when concerns overlap, one specification MUST own the requirement and the others MAY summarize the relationship.

Changes to specification paths or ownership MUST update affected agent guidance in the same change when that guidance refers to the changed specification. Depmesh governance maintenance for specification changes is owned by `./specs/behavior/files_relations.md`.

An implementation change governed by a specification MUST preserve that specification or update the specification in the same change when the intended contract changes.

## Review and maintenance

Specification changes MUST be reviewed for internal consistency, consistency with the specification index, and compatibility with existing normative statements.

New or modified specifications MUST be checked for the mandatory headings, RFC 2119 usage, stable path conventions, and corresponding Depmesh relations where applicable.

Historical implementation behavior MUST NOT override an explicit current specification. When code and specification disagree, the discrepancy MUST be reported and resolved deliberately rather than silently normalizing one to the other.
