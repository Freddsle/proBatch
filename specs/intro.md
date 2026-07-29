# List of specifications

## Goal of the document

This document lists the specification directories and specification documents in proBatch and briefly describes their purpose.

## Scope

The scope of this specification is limited to the specification index.

Detailed requirements owned by individual specifications are out of scope except for the summaries needed to keep this index useful.

## Specification directories

- `./specs/` — contains all project specifications and this index.
- `./specs/behavior/` — contains specifications for stable project behavior and artifact relationships.
- `./specs/general/` — contains cross-cutting requirements for project maintenance and development workflows.
- `./specs/meta/` — contains requirements for the structure and maintenance of specification documents.

## Specification documents

- `./specs/intro.md` — this file; indexes every specification directory and document.
- `./specs/behavior/files_relations.md` — defines the required file-relation vocabulary, directionality, discovery rules, maintenance expectations, and Depmesh behavior.
- `./specs/general/commits.md` — defines the mandatory Conventional Commits format for agent-generated commit messages.
- `./specs/general/workflows.md` — catalogs the Donna workflows required to develop, verify, document, and release the package.
- `./specs/meta/general.md` — defines general structure, style, abstraction, and maintenance requirements for all project specifications.

## Index maintenance

The specification-document list above is the complete current set. Specification-index completeness and index-maintenance requirements are owned by `./specs/meta/general.md`.

Index entries MUST describe ownership without duplicating the detailed requirements of the indexed specification.
