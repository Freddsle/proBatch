# Bootstrap Project Specifications

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "review_reference_specs"
```

Establish the specification system, file-relation model, and workflow catalog for this R/Bioconductor package.

## Review Reference Specifications

```toml donna
id = "review_reference_specs"
kind = "donna.lib.request_action"
```

1. Read the example specification index at `https://github.com/Tiendil/donna/blob/main/specs/intro.md`.
2. Read the example specification rules at `https://github.com/Tiendil/donna/blob/main/specs/meta/general.md`.
3. Inspect the editable repository structure as an R/Bioconductor package, including `DESCRIPTION`, `R/`, `tests/`, `vignettes/`, Donna configuration, and Depmesh configuration. Treat `man/` as opaque maintainer-generated output and do not access it.
4. Record the specification structure, normative style, and stable artifact families that the local documents must cover.
5. When the examples and project structure are understood, `{{ donna.lib.goto("write_core_specs") }}`.
6. If the reference material cannot be obtained or the project structure is ambiguous, `{{ donna.lib.goto("blocked") }}`.

## Write Core Specifications

```toml donna
id = "write_core_specs"
kind = "donna.lib.request_action"
```

1. Create `{{ donna.lib.path("@/specs/intro.md") }}` as the complete index of local specification directories and documents.
2. Create `{{ donna.lib.path("@/specs/meta/general.md") }}` using the reference meta-specification as a model and adapting it to this project.
3. Update `{{ donna.lib.path("@/AGENTS.md") }}` so agents must read and follow the specification hierarchy and understand the R/Bioconductor package layout.
4. Ensure both specifications contain the mandatory `Goal of the document` and `Scope` sections.
5. When these deliverables are complete, `{{ donna.lib.goto("implement_file_relations") }}`.
6. If completion is blocked, `{{ donna.lib.goto("blocked") }}`.

## Implement File Relations

```toml donna
id = "implement_file_relations"
kind = "donna.lib.request_action"
```

1. Create `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}` with detailed, normative requirements for how file relations are named, directed, discovered, maintained, and validated.
2. Update `{{ donna.lib.path("@/depmesh.toml") }}` to implement the specified relations for specifications, R sources, tests, Donna workflows, and their governing files without accessing `man/`.
3. Add project-local helper scripts under `@/bin/depemesh/` only where pattern rules cannot discover accurate relations.
4. Keep every relation bidirectional when the reverse lookup is useful, and use root-anchored artifact IDs.
5. Query representative artifacts in both directions before continuing.
6. When file relations work as specified, `{{ donna.lib.goto("document_workflows") }}`.
7. If completion is blocked, `{{ donna.lib.goto("blocked") }}`.

## Document Required Workflows

```toml donna
id = "document_workflows"
kind = "donna.lib.request_action"
```

1. Create `{{ donna.lib.path("@/specs/general/workflows.md") }}`.
2. List every Donna workflow required for maintaining this R/Bioconductor package.
3. For each workflow, define its purpose, triggers, major stages, expected verification, and completion outcome.
4. Include the current specification-bootstrap workflow in the catalog.
5. Update the specification index and any agent guidance affected by the completed workflow catalog.
6. When the catalog is complete and internally consistent, `{{ donna.lib.goto("verify_deliverables") }}`.
7. If completion is blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify Deliverables

```toml donna
id = "verify_deliverables"
kind = "donna.lib.run_script"
save_stdout_to = "verification_stdout"
save_stderr_to = "verification_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_deliverables"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

required_files=(
    "specs/intro.md"
    "specs/meta/general.md"
    "specs/behavior/files_relations.md"
    "specs/general/commits.md"
    "specs/general/workflows.md"
    "AGENTS.md"
    "depmesh.toml"
)

for file in "${required_files[@]}"; do
    test -f "$file"
done

while IFS= read -r spec; do
    grep -q '^# ' "$spec"
    grep -q '^## Goal of the document$' "$spec"
    grep -q '^## Scope$' "$spec"
done < <(find specs -type f -name '*.md' -print)

grep -q 'specs/intro.md' AGENTS.md
grep -q 'specs/general/commits.md' AGENTS.md
grep -q 'R/' AGENTS.md
grep -q 'Bioconductor' AGENTS.md

depmesh -p llm relations
depmesh -p llm dependencies --relation governed_by @/depmesh.toml
depmesh -p llm dependencies --relation governs @/specs/behavior/files_relations.md
depmesh -p llm dependencies --relation tested_by @/R/normalize.R
depmesh -p llm dependencies --relation tests @/tests/testthat/test-normalize.R
if grep -q '@/man/' depmesh.toml; then
    echo "depmesh.toml must not reference man/" >&2
    exit 1
fi

donna -p llm validate @/workflows/bootstrap-project-specifications.donna.md
```

## Repair Deliverables

```toml donna
id = "repair_deliverables"
kind = "donna.lib.request_action"
```

The verification step failed.

Standard output:

```text
{{ donna.lib.task_variable("verification_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("verification_stderr") }}
```

Diagnose and repair the deliverables, then `{{ donna.lib.goto("verify_deliverables") }}`.

If the failure cannot be repaired without developer input, `{{ donna.lib.goto("blocked") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The specification system, file relations, and Donna workflow catalog are complete and verified. Report the workflow execution result, files created or changed, relation queries performed, and any intentionally deferred workflows.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The workflow cannot continue without developer input. Report the completed work, exact blocker, and required decision.
