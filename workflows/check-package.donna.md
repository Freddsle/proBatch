# Check R and Bioconductor Package

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Run the authoritative agent-safe portion of the local R/Bioconductor package quality gate, compose the focused documentation, test, and vignette workflows when their triggers apply, and hand every package-wide or generated-document check to the maintainer.

This workflow MUST NOT access `man/`, generate documentation, modify generated `NAMESPACE`, run a package-wide build, package check, or Bioconductor check, install or update dependencies or the Pixi environment, use the network, stage files, create commits, publish artifacts, or reset the Donna session.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/meta/general.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}` completely.
2. Read the built-in Donna workflow and Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations separately for this workflow, its catalog, `{{ donna.lib.path("@/AGENTS.md") }}`, and every non-generated package artifact already known to be in the initiating scope.
3. Inspect the worktree without modifying the Git index. Exclude `man/` from every path query and diff. Preserve unrelated user changes and distinguish the initiating package-check scope from pre-existing changes or failures.
4. Confirm this action belongs to the current Donna task and no unrelated action request is being displaced. Do not reset the session.
5. Confirm this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna operation was executing.
6. Confirm that agent-safe metadata, source, test, and isolated-vignette checks are authorized. Package installation, documentation generation, package-wide building or checking, Bioconductor checking, dependency installation, network access, staging, commits, publication, and submission remain unauthorized.
7. Do not read, enumerate, compare, hash, validate, or otherwise inspect any path under `man/`. Do not modify or regenerate `NAMESPACE`. Leave every command that may inspect generated manuals to the maintainer.
8. If the governed scope, current diff, and agent/maintainer boundary are understood, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
9. If specifications conflict, required context is unavailable, or the requested scope requires unauthorized work, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "verify_package_environment"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/check-package.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"

artifact_entry_count="$(
    rg -Fxc -- '- Artifact: `./workflows/check-package.donna.md`' \
        specs/general/workflows.md
)"
if (( artifact_entry_count != 1 )); then
    printf 'Workflow catalog must contain exactly one check-package artifact entry; found %d.\n' \
        "$artifact_entry_count" >&2
    exit 1
fi

catalog_state="$(
    awk '
        /^### Check R and Bioconductor package$/ {
            in_section = 1
            next
        }
        in_section && /^### / {
            exit
        }
        in_section && /^- State: / {
            sub(/^- State: /, "")
            print
            exit
        }
    ' specs/general/workflows.md
)"
case "$catalog_state" in
    "implementation in progress"|"implemented") ;;
    *)
        printf 'Unexpected check-package catalog state: %s\n' \
            "$catalog_state" >&2
        exit 1
        ;;
esac

if ! rg -Fq -- '@/workflows/check-package.donna.md' AGENTS.md; then
    printf 'AGENTS.md must route package-check triggers to this workflow.\n' >&2
    exit 1
fi

for child_workflow in \
    workflows/review-documentation-sources.donna.md \
    workflows/run-tests.donna.md \
    workflows/build-vignettes.donna.md
do
    if [[ ! -f "$child_workflow" ]]; then
        printf 'Required child workflow is missing: %s\n' \
            "$child_workflow" >&2
        exit 1
    fi
done

relation_artifacts() {
    local relation="$1"
    local artifact="$2"
    local output
    output="$(depmesh -p llm dependencies --relation "$relation" "$artifact")"
    if printf '%s\n' "$output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for %s on %s:\n%s\n' \
            "$relation" "$artifact" "$output" >&2
        return 1
    fi
    printf '%s\n' "$output" |
        sed -n 's/^- \(@\/[^[:space:]]\+\)$/\1/p' |
        LC_ALL=C sort
}

actual_governance="$(relation_artifacts governed_by "$workflow_id")"
expected_governance="$(
    printf '%s\n' \
        @/specs/behavior/files_relations.md \
        @/specs/general/workflows.md |
        LC_ALL=C sort
)"
if [[ "$actual_governance" != "$expected_governance" ]]; then
    printf 'Unexpected workflow governance.\nExpected:\n%s\nActual:\n%s\n' \
        "$expected_governance" "$actual_governance" >&2
    exit 1
fi

for governing_spec in \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
do
    if ! relation_artifacts governs "$governing_spec" |
        rg -Fqx -- "$workflow_id"; then
        printf 'Reverse governance does not contain %s for %s.\n' \
            "$workflow_id" "$governing_spec" >&2
        exit 1
    fi
done

if ! rg -Fqx -- '^workflows$' .Rbuildignore; then
    printf 'Project workflows must remain excluded from the source package.\n' >&2
    exit 1
fi

safe_diff_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks diff --check -- "${safe_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_diff_paths[@]}"

for text_path in \
    "$workflow_path" \
    AGENTS.md \
    specs/general/workflows.md
do
    if rg -n '[[:blank:]]+$' "$text_path"; then
        printf '%s contains trailing whitespace.\n' "$text_path" >&2
        exit 1
    fi
    if LC_ALL=C rg -n $'\r$' "$text_path"; then
        printf '%s contains CRLF line endings.\n' "$text_path" >&2
        exit 1
    fi
    if [[ -n "$(tail -c 1 -- "$text_path")" ]]; then
        printf '%s must end with a newline.\n' "$text_path" >&2
        exit 1
    fi
done

printf 'Workflow source, lifecycle state, routing, child workflows, governance, package exclusion, and non-generated diff hygiene are synchronized.\n'
```

## Repair Workflow Artifact

```toml donna
id = "repair_workflow_artifact"
kind = "donna.lib.request_action"
```

Workflow-artifact verification failed.

Standard output:

```text
{{ donna.lib.task_variable("workflow_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("workflow_stderr") }}
```

1. Query all configured relations separately for every artifact to be edited.
2. Repair only this workflow, its catalog entry, affected agent guidance, its scalable governance relation, or an in-scope formatting discrepancy.
3. Preserve unrelated changes and the Git index. Do not add a workflow-specific Depmesh rule when the scalable workflow rule applies.
4. Keep the catalog state at `implementation in progress` until the first representative execution finishes.
5. Validate every control-flow change while Donna is awaiting this action.
6. After repair, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
7. If repair requires a new contract or developer decision, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Package Environment

```toml donna
id = "verify_package_environment"
kind = "donna.lib.run_script"
save_stdout_to = "environment_stdout"
save_stderr_to = "environment_stderr"
goto_on_success = "review_package_scope"
goto_on_failure = "diagnose_environment"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

test -f DESCRIPTION
test -f .Rbuildignore
test -f pixi.toml
test -f pixi.lock
test -f tests/testthat.R
test -d R
test -d tests/testthat
command -v pixi >/dev/null 2>&1
command -v depmesh >/dev/null 2>&1
command -v donna >/dev/null 2>&1

state_dir=".session/donna/check-package"
if [[ -L "$state_dir" ]]; then
    printf 'Package-check state directory must not be a symlink.\n' >&2
    exit 1
fi
mkdir -p "$state_dir"

safe_index_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks rev-parse --verify HEAD >"$state_dir/head"
git --no-optional-locks diff --cached --binary -- \
    "${safe_index_paths[@]}" |
    git hash-object --stdin >"$state_dir/index-hash"

unexpected_artifacts="$(
    find . -maxdepth 1 \
        \( \
            -type d -name '*.Rcheck' \
            -o -type f -name '*.tar.gz' \
            -o -type f -name 'Rplots.pdf' \
        \) -print
)"
if [[ -n "$unexpected_artifacts" ]]; then
    printf 'Unexpected pre-check build or runtime artifacts exist:\n' >&2
    printf '%s\n' "$unexpected_artifacts" >&2
    exit 1
fi

pixi run --as-is Rscript --vanilla - <<'RSCRIPT'
description <- read.dcf("DESCRIPTION")
if (nrow(description) != 1L) {
    stop("DESCRIPTION must contain exactly one DCF record.", call. = FALSE)
}
required <- c(
    "Package", "Title", "Version", "Authors@R", "Description",
    "License", "Depends", "Imports", "Suggests", "Encoding"
)
missing <- setdiff(required, colnames(description))
if (length(missing) > 0L) {
    stop(
        sprintf(
            "DESCRIPTION is missing required fields: %s",
            paste(missing, collapse = ", ")
        ),
        call. = FALSE
    )
}
if (!identical(unname(description[1L, "Package"]), "proBatch")) {
    stop("DESCRIPTION Package must remain proBatch.", call. = FALSE)
}
cat(R.version.string, "\n", sep = "")
cat(
    "Parsed DESCRIPTION for proBatch ",
    unname(description[1L, "Version"]),
    ".\n",
    sep = ""
)
RSCRIPT

printf 'Package-check environment and protected Git baseline are ready.\n'
```

## Diagnose Environment

```toml donna
id = "diagnose_environment"
kind = "donna.lib.request_action"
```

Package-environment verification failed.

Standard output:

```text
{{ donna.lib.task_variable("environment_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("environment_stderr") }}
```

1. Determine whether the failure is a missing project artifact, malformed editable metadata, an unexpected build/runtime artifact, a missing existing local tool, or a verifier defect.
2. Query all configured relations before editing any governed artifact. Preserve unrelated changes and the Git index.
3. Repair only an authorized editable source or this workflow, then `{{ donna.lib.goto("verify_package_environment") }}`.
4. Do not install or update dependencies, use the network, inspect generated manuals, generate documentation, or run a package-wide command.
5. If the environment cannot be made ready within the authorized local scope, `{{ donna.lib.goto("record_blocker") }}`.

## Review Package Scope

```toml donna
id = "review_package_scope"
kind = "donna.lib.request_action"
```

Environment verification output:

```text
{{ donna.lib.task_variable("environment_stdout") }}
```

1. Inspect `DESCRIPTION`, `pixi.toml`, `pixi.lock`, `.Rbuildignore`, and the initiating non-generated diff without accessing `man/` or diff-reviewing generated `NAMESPACE`.
2. Identify changes to dependency fields, supported R, build metadata, R or Roxygen2 sources, exports or registrations at their editable source, tests, examples, package data, vignette sources/assets, and release intent.
3. Query all configured Depmesh relations separately for every artifact that may need review or repair. Do not interpret an empty result as proof of no conceptual dependency.
4. Confirm that the complete testthat suite is required for this quality gate even when no focused test is selected.
5. Run documentation-source review when Roxygen2 comments, public signatures, exports, registrations, or documented datasets changed, or when release scope requires it.
6. Run the vignette workflow when vignettes, their assets, APIs/data they execute, documentation dependencies, or release scope changed.
7. If documentation-source review is required, `{{ donna.lib.goto("run_documentation_review") }}`.
8. If documentation-source review is not required, `{{ donna.lib.goto("run_test_suite") }}`.
9. If the scope is ambiguous or conflicts with project requirements, `{{ donna.lib.goto("record_blocker") }}`.

## Run Documentation Review

```toml donna
id = "run_documentation_review"
kind = "donna.lib.request_action"
```

1. Run `{{ donna.lib.path("@/workflows/review-documentation-sources.donna.md") }}` as a child workflow and complete it before returning to this request.
2. Keep this parent action request pending while the child workflow is active.
3. Treat the child as successful here only when its source-level verification passed. Its normal handoff of generated documentation to the maintainer does not invalidate that result; a handoff without completed source gates or a blocked finish is not successful.
4. If the child workflow completes its source-level verification successfully, `{{ donna.lib.goto("run_test_suite") }}`.
5. If it identifies an in-scope repair that this parent must reconsider, `{{ donna.lib.goto("repair_editable_sources") }}`.
6. If it is blocked by generated output, environment, authority, or an unresolved contract, `{{ donna.lib.goto("record_blocker") }}`.

## Run Complete Test Suite

```toml donna
id = "run_test_suite"
kind = "donna.lib.request_action"
```

1. Run `{{ donna.lib.path("@/workflows/run-tests.donna.md") }}` as a child workflow and complete it before returning to this request.
2. Tell the child that this package quality gate requires the complete testthat suite; focused tests remain useful when the initiating scope identifies them.
3. Keep this parent action request pending while the child workflow is active.
4. If the child workflow finishes with affected focused tests and the complete suite passing, `{{ donna.lib.goto("decide_vignette_build") }}`.
5. If it finishes blocked or reports a test/environment limitation that prevents the required complete suite, `{{ donna.lib.goto("record_blocker") }}`.

## Decide Vignette Build

```toml donna
id = "decide_vignette_build"
kind = "donna.lib.request_action"
```

1. Recheck the initiating scope and any repairs made by the documentation or test child workflow.
2. If vignettes, vignette assets, APIs/data executed by them, documentation dependencies, or release scope are affected, `{{ donna.lib.goto("run_vignette_build") }}`.
3. If none of those triggers applies, record why the specialized vignette workflow is not required and `{{ donna.lib.goto("verify_editable_sources") }}`.
4. If the trigger decision cannot be made from the available evidence, `{{ donna.lib.goto("record_blocker") }}`.

## Run Vignette Build

```toml donna
id = "run_vignette_build"
kind = "donna.lib.request_action"
```

1. Run `{{ donna.lib.path("@/workflows/build-vignettes.donna.md") }}` as a child workflow and complete it before returning to this request.
2. Keep this parent action request pending while the child workflow is active.
3. If every declared vignette builds successfully and the child confirms isolated runtime-state hygiene, `{{ donna.lib.goto("verify_editable_sources") }}`.
4. If the child workflow finishes blocked or requires dependency installation, network access, generated documentation, or work outside the authorized scope, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Editable Sources

```toml donna
id = "verify_editable_sources"
kind = "donna.lib.run_script"
save_stdout_to = "source_stdout"
save_stderr_to = "source_stderr"
goto_on_success = "review_agent_results"
goto_on_failure = "repair_editable_sources"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

safe_diff_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks diff --check -- "${safe_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_diff_paths[@]}"

for required_exclusion in \
    '^AGENTS\.md$' \
    '^bin$' \
    '^depmesh\.toml$' \
    '^donna\.toml$' \
    '^specs$' \
    '^workflows$' \
    '^\.session$'
do
    if ! rg -Fqx -- "$required_exclusion" .Rbuildignore; then
        printf 'Missing project-only .Rbuildignore exclusion: %s\n' \
            "$required_exclusion" >&2
        exit 1
    fi
done

pixi run --as-is Rscript --vanilla - <<'RSCRIPT'
description <- read.dcf("DESCRIPTION")

field_packages <- function(field) {
    if (!field %in% colnames(description)) {
        return(character())
    }
    entries <- trimws(strsplit(description[1L, field], ",", fixed = TRUE)[[1L]])
    entries <- sub("[[:space:]]*\\(.*\\)[[:space:]]*$", "", entries)
    entries[nzchar(entries)]
}

imports <- field_packages("Imports")
suggests <- field_packages("Suggests")
depends <- field_packages("Depends")

for (field_name in c("Imports", "Suggests", "Depends")) {
    packages <- field_packages(field_name)
    duplicated_packages <- unique(packages[duplicated(packages)])
    if (length(duplicated_packages) > 0L) {
        stop(
            sprintf(
                "%s contains duplicate dependencies: %s",
                field_name,
                paste(duplicated_packages, collapse = ", ")
            ),
            call. = FALSE
        )
    }
}

overlap <- intersect(imports, suggests)
if (length(overlap) > 0L) {
    stop(
        sprintf(
            "Dependencies occur in both Imports and Suggests: %s",
            paste(overlap, collapse = ", ")
        ),
        call. = FALSE
    )
}

if (!"R" %in% depends) {
    stop("DESCRIPTION Depends must declare the supported R version.", call. = FALSE)
}

if ("VignetteBuilder" %in% colnames(description)) {
    builders <- field_packages("VignetteBuilder")
    undeclared <- setdiff(builders, c(imports, suggests))
    if (length(undeclared) > 0L) {
        stop(
            sprintf(
                "VignetteBuilder packages are undeclared: %s",
                paste(undeclared, collapse = ", ")
            ),
            call. = FALSE
        )
    }
}

source_paths <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
if (length(source_paths) == 0L) {
    stop("No package R sources were found.", call. = FALSE)
}

example_paths <- sort(
    list.files("inst/examples", pattern = "[.]R$", full.names = TRUE)
)
for (source_path in c(source_paths, example_paths)) {
    parse(file = source_path, keep.source = TRUE)
}

cat(
    sprintf(
        "Parsed DESCRIPTION, %d package R sources, and %d installed R examples.\n",
        length(source_paths),
        length(example_paths)
    )
)
cat(
    sprintf(
        "Dependency fields contain %d Imports and %d Suggests without duplicates or overlap.\n",
        length(imports),
        length(suggests)
    )
)
RSCRIPT

relations_output="$(depmesh -p llm relations)"
for relation in governed_by governs tested_by tests; do
    if ! printf '%s\n' "$relations_output" |
        rg -Fqx -- "## ${relation}"; then
        printf 'Depmesh relation is missing: %s\n' "$relation" >&2
        exit 1
    fi
done

workflow_governance="$(
    depmesh -p llm dependencies \
        --relation governed_by \
        @/workflows/check-package.donna.md
)"
for governing_spec in \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
do
    if ! printf '%s\n' "$workflow_governance" |
        rg -Fqx -- "- ${governing_spec}"; then
        printf 'Workflow governance is missing: %s\n' \
            "$governing_spec" >&2
        exit 1
    fi
done

unexpected_artifacts="$(
    find . -maxdepth 1 \
        \( \
            -type d -name '*.Rcheck' \
            -o -type f -name '*.tar.gz' \
            -o -type f -name 'Rplots.pdf' \
        \) -print
)"
if [[ -n "$unexpected_artifacts" ]]; then
    printf 'Agent-safe checks left unexpected build or runtime artifacts:\n' >&2
    printf '%s\n' "$unexpected_artifacts" >&2
    exit 1
fi

printf 'Editable metadata, R syntax, dependency fields, project exclusions, Depmesh availability, and runtime-artifact hygiene pass agent-safe checks.\n'
```

## Repair Editable Sources

```toml donna
id = "repair_editable_sources"
kind = "donna.lib.request_action"
```

Agent-safe editable-source verification or a child review found a repairable issue.

Standard output:

```text
{{ donna.lib.task_variable("source_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("source_stderr") }}
```

1. Classify the finding as editable metadata, dependency declarations, R or Roxygen2 source, tests, examples, vignette source/assets, runtime artifacts, environment, or a verifier defect.
2. Query all configured relations separately for every artifact before editing it. Preserve unrelated changes and the Git index.
3. Repair only authorized non-generated sources. Add or update focused testthat coverage for executable R changes and preserve public API documentation conventions at their Roxygen2 source.
4. Do not inspect or repair `man/`, modify generated `NAMESPACE`, generate documentation, install dependencies, use the network, or run a package-wide command.
5. After repairing any editable package metadata, dependency declaration, R or Roxygen2 source, test, example, package-data artifact, or vignette source/asset, `{{ donna.lib.goto("review_package_scope") }}` so every applicable child workflow and deterministic gate runs again.
6. If the finding is limited to this workflow's deterministic verifier, repair and validate it while Donna is awaiting this action, then `{{ donna.lib.goto("verify_workflow_artifact") }}`.
7. If the issue cannot be repaired within the initiating scope, `{{ donna.lib.goto("record_blocker") }}`.

## Review Agent Results

```toml donna
id = "review_agent_results"
kind = "donna.lib.request_action"
```

Agent-safe verification output:

```text
{{ donna.lib.task_variable("source_stdout") }}
```

1. Confirm the complete testthat child workflow passed and report its focused-test selection, skips, and environment limitations.
2. Confirm each documentation or vignette child workflow required by the reviewed scope completed successfully, or record the evidence that its trigger did not apply.
3. Review the editable metadata and dependency decisions against `DESCRIPTION`, Pixi metadata, the initiating scope, and project requirements.
4. Confirm the deterministic source check parsed all editable package R sources and installed R examples, and that no package build, package check, Bioconductor check, documentation generation, dependency installation, network access, or generated-manual inspection occurred.
5. Confirm no unexpected build or runtime artifact remains.
6. If every agent-owned gate passes, `{{ donna.lib.goto("prepare_maintainer_handoff") }}`.
7. If an in-scope discrepancy remains, `{{ donna.lib.goto("repair_editable_sources") }}`.
8. If acceptance requires a product decision, unavailable environment, generated-output inspection, or broader authority, `{{ donna.lib.goto("record_blocker") }}`.

## Prepare Maintainer Handoff

```toml donna
id = "prepare_maintainer_handoff"
kind = "donna.lib.request_action"
```

The agent MUST NOT perform the maintainer-owned operations in this request.

1. Prepare an exact handoff asking the maintainer to regenerate Roxygen2 documentation and `NAMESPACE` manually with devtools in the supported local environment.
2. Ask the maintainer to run the package-wide R checks and Bioconductor checks appropriate to the initiating scope. These checks may inspect generated manuals and therefore remain maintainer-owned.
3. Ask for the R/Bioconductor environment, commands, exit outcomes, and all errors, warnings, or notes, classified as editable-source, generated-manual/namespace, environment, or policy findings.
4. State that the agent will not inspect, diff-review, lint, validate, or repair `man/`; a generated-manual problem remains entirely with the maintainer.
5. If the maintainer has already supplied complete passing confirmation for generation and every required package-wide check, `{{ donna.lib.goto("verify_runtime_state") }}`.
6. If no confirmation is available yet, `{{ donna.lib.goto("finish_handoff") }}` and report the handoff without claiming the authoritative package gate passed.
7. If the maintainer reports a non-generated editable-source failure within scope, `{{ donna.lib.goto("repair_maintainer_failure") }}`.
8. If the maintainer reports only a generated-manual or generated-namespace problem, `{{ donna.lib.goto("finish_handoff") }}` and defer it without inspection or repair.
9. If the required check environment or a product decision is unavailable, `{{ donna.lib.goto("record_blocker") }}`.

## Repair Maintainer-Reported Failure

```toml donna
id = "repair_maintainer_failure"
kind = "donna.lib.request_action"
```

1. Use only the maintainer's report; do not inspect generated manuals or generated `NAMESPACE`.
2. Classify the failure and query all configured relations for every editable artifact that may need repair.
3. If the issue is an in-scope non-generated source, metadata, test, example, or vignette defect, repair it and `{{ donna.lib.goto("review_package_scope") }}` so every applicable agent-owned gate repeats.
4. If the issue concerns `man/`, generated `NAMESPACE`, or documentation generation itself, leave it entirely to the maintainer and `{{ donna.lib.goto("finish_handoff") }}`.
5. If repair requires scope expansion, dependency installation, network access, or a new product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Runtime State

```toml donna
id = "verify_runtime_state"
kind = "donna.lib.run_script"
save_stdout_to = "runtime_stdout"
save_stderr_to = "runtime_stderr"
goto_on_success = "final_review"
goto_on_failure = "repair_runtime_state"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

state_dir=".session/donna/check-package"
test -f "$state_dir/head"
test -f "$state_dir/index-hash"

expected_head="$(<"$state_dir/head")"
actual_head="$(git --no-optional-locks rev-parse --verify HEAD)"
if [[ "$actual_head" != "$expected_head" ]]; then
    printf 'HEAD changed during the package-check workflow.\n' >&2
    exit 1
fi

safe_index_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
expected_index_hash="$(<"$state_dir/index-hash")"
actual_index_hash="$(
    git --no-optional-locks diff --cached --binary -- \
        "${safe_index_paths[@]}" |
        git hash-object --stdin
)"
if [[ "$actual_index_hash" != "$expected_index_hash" ]]; then
    printf 'The protected Git index changed during the package-check workflow.\n' >&2
    exit 1
fi

git --no-optional-locks diff --check -- "${safe_index_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_index_paths[@]}"

unexpected_artifacts="$(
    find . -maxdepth 1 \
        \( \
            -type d -name '*.Rcheck' \
            -o -type f -name '*.tar.gz' \
            -o -type f -name 'Rplots.pdf' \
        \) -print
)"
if [[ -n "$unexpected_artifacts" ]]; then
    printf 'Maintainer checks left unexpected package artifacts:\n' >&2
    printf '%s\n' "$unexpected_artifacts" >&2
    exit 1
fi

printf 'HEAD, protected index state, non-generated diff hygiene, and package-artifact hygiene remain valid after maintainer confirmation.\n'
```

## Repair Runtime State

```toml donna
id = "repair_runtime_state"
kind = "donna.lib.request_action"
```

Runtime-state verification failed after maintainer confirmation.

Standard output:

```text
{{ donna.lib.task_variable("runtime_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("runtime_stderr") }}
```

1. Do not discard, reset, stage, or overwrite any user or maintainer change.
2. If HEAD or the protected Git index changed, preserve the change and `{{ donna.lib.goto("record_blocker") }}` so the maintainer can decide how to proceed.
3. If a new top-level check/build artifact can be safely removed only with explicit maintainer approval, `{{ donna.lib.goto("record_blocker") }}`.
4. If the verifier alone is defective, query this workflow's relations, repair and validate it, then `{{ donna.lib.goto("verify_workflow_artifact") }}`.
5. If the runtime state again satisfies the captured baseline without agent mutation, `{{ donna.lib.goto("verify_runtime_state") }}`.

## Final Review

```toml donna
id = "final_review"
kind = "donna.lib.request_action"
```

Runtime verification output:

```text
{{ donna.lib.task_variable("runtime_stdout") }}
```

1. Confirm all required agent-owned child workflows and deterministic gates passed after the final editable-source repair.
2. Confirm the maintainer explicitly reported successful documentation generation, package-wide R checks, and required Bioconductor checks, including any accepted warnings or notes.
3. Do not inspect generated manuals or generated `NAMESPACE`; rely only on the maintainer's report for those artifacts.
4. While Donna is awaiting this action and no workflow operation is executing, run `donna -p llm validate @/workflows/check-package.donna.md`, render it in view mode, and run `donna -p llm list`; inspect the transitions and listed summary.
5. Confirm the catalog, AGENTS routing, scalable Depmesh governance, child composition, and workflow artifact form one coherent implementation.
6. If the authoritative local package quality gate is confirmed, `{{ donna.lib.goto("finish_verified") }}`.
7. If a repairable non-generated discrepancy remains, `{{ donna.lib.goto("repair_editable_sources") }}`.
8. If a generated-output, environment, contract, or authority blocker remains, `{{ donna.lib.goto("record_blocker") }}`.

## Finish Verified

```toml donna
id = "finish_verified"
kind = "donna.lib.finish"
```

The agent-owned source, metadata, test, and applicable vignette/documentation-source gates passed, and the maintainer confirmed documentation generation plus package-wide R and Bioconductor checks. Report the initiating scope, child workflow outcomes, deterministic evidence, maintainer commands and results, accepted warnings or notes, and final runtime-state verification without inspecting generated manuals.

## Finish With Maintainer Handoff

```toml donna
id = "finish_handoff"
kind = "donna.lib.finish"
```

The agent-owned gates passed and the package-wide generation and checking work has been explicitly handed to the maintainer. Report completed child workflows and deterministic checks, the exact maintainer-owned commands or outcomes still required, and any source-level findings. Do not claim the authoritative package quality gate passed until the maintainer confirms those operations.

## Record Blocker

```toml donna
id = "record_blocker"
kind = "donna.lib.request_action"
```

1. Summarize the completed gates, exact blocker, affected editable artifacts, preserved user changes, and the developer or environment decision required.
2. If no remaining workflow resolves a repository finding, record it in `{{ donna.lib.path("@/workflows/reports/package-check-open-questions.md") }}`. Query Depmesh before creating or editing the report, keep it non-executable, and do not add it to the workflow catalog.
3. Do not record transient missing maintainer confirmation as an open question; use the explicit handoff finish for that ordinary outcome.
4. Do not inspect or repair generated manuals, modify generated `NAMESPACE`, install dependencies, use the network, stage files, create commits, publish, submit, or reset the Donna session.
5. After recording any required persistent finding, `{{ donna.lib.goto("blocked") }}`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The package-check workflow cannot continue within the authorized agent-safe scope. Report completed gates, the exact source, environment, generated-output, contract, or authority blocker, any persistent open-question entry, affected editable artifacts, and the decision required. Do not claim package verification, inspect generated manuals, modify generated `NAMESPACE`, stage files, create commits, publish, submit, or reset the Donna session.
