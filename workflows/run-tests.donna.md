# Run Focused Package Tests

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Select focused testthat tests from the initiating scope and Depmesh relations, run them when useful, then run the complete testthat suite in a separate clean R process and guide only authorized repairs.

This workflow MUST NOT access `man/`, use the network, install or update the Pixi environment, run package-wide checks, generate documentation, stage files, create commits, or reset the Donna session.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}` completely.
2. Read the built-in Donna and Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations separately for every R source, test, fixture, or package-data artifact already known to be in the initiating scope.
3. Inspect the current worktree without modifying the Git index. Preserve unrelated user changes and distinguish the initiating scope from pre-existing failures.
4. Confirm that this action belongs to the current Donna task and that no unrelated action request is being displaced. Do not reset the session.
5. Confirm that this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna operation was executing.
6. Confirm that test execution is authorized but package installation, dependency installation, network access, package checks, documentation generation, staging, and commits are not.
7. Do not read, enumerate, compare, or otherwise inspect any path under `man/`.
8. If the governed scope, current diff, and test-only boundary are understood, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
9. If required context is unavailable or the requested testing scope conflicts with a specification, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "verify_test_environment"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/run-tests.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"

artifact_entry_count="$(
    rg -Fxc -- '- Artifact: `./workflows/run-tests.donna.md`' \
        specs/general/workflows.md
)"
if (( artifact_entry_count != 1 )); then
    printf 'Workflow catalog must contain exactly one run-tests artifact entry; found %d.\n' \
        "$artifact_entry_count" >&2
    exit 1
fi

catalog_state="$(
    awk '
        /^### Run focused package tests$/ {
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
        printf 'Unexpected run-tests catalog state: %s\n' "$catalog_state" >&2
        exit 1
        ;;
esac

if ! rg -Fq -- '@/workflows/run-tests.donna.md' AGENTS.md; then
    printf 'AGENTS.md must route testthat regression triggers to this workflow.\n' >&2
    exit 1
fi

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

git diff --check -- . ':(top,exclude)man/**'
git diff --cached --check -- . ':(top,exclude)man/**'

if rg -n '[[:blank:]]+$' "$workflow_path"; then
    printf 'Workflow contains trailing whitespace.\n' >&2
    exit 1
fi
if LC_ALL=C rg -n $'\r$' "$workflow_path"; then
    printf 'Workflow contains CRLF line endings.\n' >&2
    exit 1
fi
if [[ -n "$(tail -c 1 "$workflow_path")" ]]; then
    printf 'Workflow must end with a newline.\n' >&2
    exit 1
fi

printf 'Workflow source, catalog state, agent guidance, package exclusion, governance, and non-man diff hygiene are synchronized.\n'
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
3. Preserve unrelated changes and the Git index. Do not weaken the generated-manual exclusion or add a workflow-specific Depmesh rule when the scalable workflow rule applies.
4. Keep the catalog state at `implementation in progress` until the first representative execution finishes.
5. After repair, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
6. If repair requires a new contract or developer decision, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Test Environment

```toml donna
id = "verify_test_environment"
kind = "donna.lib.run_script"
save_stdout_to = "environment_stdout"
save_stderr_to = "environment_stderr"
goto_on_success = "discover_test_plan"
goto_on_failure = "diagnose_environment"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

test -f DESCRIPTION
test -f pixi.toml
test -f pixi.lock
test -f tests/testthat.R
test -d tests/testthat
command -v pixi >/dev/null 2>&1

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected pre-test graphics artifact: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

assert_graphics_artifacts_not_ignored() {
    local ignore_output ignore_status
    set +e
    ignore_output="$(
        git check-ignore --no-index -v -- \
            Rplots.pdf tests/testthat/Rplots.pdf 2>&1
    )"
    ignore_status=$?
    set -e
    case "$ignore_status" in
        0)
            printf 'Test graphics output is hidden by an ignore rule:\n%s\n' \
                "$ignore_output" >&2
            return 1
            ;;
        1)
            return 0
            ;;
        *)
            printf 'Unable to verify test-graphics ignore state:\n%s\n' \
                "$ignore_output" >&2
            return "$ignore_status"
            ;;
    esac
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

assert_absent_graphics_artifacts
assert_graphics_artifacts_not_ignored

plan_dir=".session/donna/run-tests"
mkdir -p "$plan_dir"
sha256sum pixi.toml pixi.lock > "$plan_dir/environment.sha256"
hash_non_man_worktree > "$plan_dir/worktree.sha256"

pixi run --as-is Rscript --vanilla -e '
options(mc.cores = 1L)
description <- read.dcf("DESCRIPTION")
if (!identical(unname(description[1, "Config/testthat/edition"]), "3")) {
    stop("Config/testthat/edition must be 3")
}
required <- c("testthat", "pkgload", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    stop("Missing test dependencies in the existing Pixi environment: ",
         paste(missing, collapse = ", "))
}
if (utils::packageVersion("testthat") < "3.0.0") {
    stop("testthat 3.0.0 or newer is required")
}
cat("R_VERSION=", R.version.string, "\n", sep = "")
cat("TESTTHAT_VERSION=", as.character(utils::packageVersion("testthat")), "\n", sep = "")
cat("TESTTHAT_EDITION=", description[1, "Config/testthat/edition"], "\n", sep = "")
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
if (!inherits(BiocParallel::bpparam(), "SerialParam")) {
    stop("The default BiocParallel backend is not SerialParam")
}
cat("BIOCPARALLEL_BACKEND=SerialParam\n")
cat("BIOCPARALLEL_PREREQUISITE=available in existing Pixi environment\n")
'

sha256sum --check "$plan_dir/environment.sha256"
assert_absent_graphics_artifacts
if [[ "$(hash_non_man_worktree)" != "$(cat "$plan_dir/worktree.sha256")" ]]; then
    printf 'The Git-visible non-man worktree differs from its environment-stage snapshot.\n' >&2
    exit 1
fi
printf 'The existing Pixi environment is ready for non-installing, clean-process testthat execution; pixi.toml and pixi.lock are unchanged.\n'
```

## Diagnose Environment

```toml donna
id = "diagnose_environment"
kind = "donna.lib.request_action"
```

Test-environment verification failed.

Standard output:

```text
{{ donna.lib.task_variable("environment_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("environment_stderr") }}
```

1. Determine whether the failure is a missing existing executable or package, stale or mutated lock metadata, a pre-existing graphics artifact, or another local limitation.
2. Do not install packages, hydrate an environment, update a lock, or use the network unless the initiating request separately authorizes that action.
3. If the existing local environment can be repaired without network access or out-of-scope changes, make the focused repair and `{{ donna.lib.goto("verify_test_environment") }}`.
4. If provisioning, generated documentation, package-wide checking, or developer authority is required, `{{ donna.lib.goto("record_blocker") }}`.

## Discover Test Plan

```toml donna
id = "discover_test_plan"
kind = "donna.lib.run_script"
save_stdout_to = "discovery_stdout"
save_stderr_to = "discovery_stderr"
goto_on_success = "review_test_plan"
goto_on_failure = "repair_test_plan_discovery"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/run-tests"
scope_file="$plan_dir/scope.txt"
focus_file="$plan_dir/focused-tests.txt"
plan_file="$plan_dir/plan.md"
changed_file="$plan_dir/changed-relative.tmp"
focus_tmp="$plan_dir/focused-tests.tmp"

mkdir -p "$plan_dir"
: > "$changed_file"
: > "$focus_tmp"

{
    git diff --name-only --diff-filter=ACDMRTUXB -- R tests data inst
    git diff --cached --name-only --diff-filter=ACDMRTUXB -- R tests data inst
    git ls-files --others --exclude-standard -- R tests data inst
} |
    sed '/^$/d' |
    LC_ALL=C sort -u > "$changed_file"

sed 's#^#@/#' "$changed_file" > "$scope_file"

{
    printf '# Run-tests scope plan\n\n'
    printf '## Scope rationale\n\n'
    if [[ -s "$scope_file" ]]; then
        printf 'Candidate scope was derived from staged, unstaged, and untracked paths under `R/`, `tests/`, `data/`, and `inst/`. Agent review must remove unrelated or pre-existing entries before this file becomes the authorized repair boundary.\n'
    else
        printf 'No changed R, test, fixture, installed-resource, or package-data artifact was discoverable from the worktree. The complete suite remains required; focused tests may remain empty unless the initiating request supplies additional scope.\n'
    fi
    printf '\n## Dependency queries\n'
} > "$plan_file"

broad_scope=0
deleted_test_scope=0
while IFS= read -r artifact_id; do
    [[ -n "$artifact_id" ]] || continue

    query_output="$(depmesh -p llm dependencies "$artifact_id")"
    if printf '%s\n' "$query_output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for %s:\n%s\n' \
            "$artifact_id" "$query_output" >&2
        exit 1
    fi

    {
        printf '\n### `%s`\n\n' "$artifact_id"
        printf '```text\n%s\n```\n' "$query_output"
    } >> "$plan_file"

    case "$artifact_id" in
        @/R/*.R)
            printf '%s\n' "$query_output" |
                sed -n 's/^- \(@\/tests\/testthat\/test-[A-Za-z0-9_-]\+\.R\)$/\1/p' \
                >> "$focus_tmp"
            ;;
        @/tests/testthat/test-*.R)
            test_path="${artifact_id#@/}"
            if [[ -f "$test_path" ]]; then
                printf '%s\n' "$artifact_id" >> "$focus_tmp"
            else
                deleted_test_scope=1
                broad_scope=1
            fi
            ;;
        @/tests/testthat/helper-*|@/tests/testthat.R|@/tests/spelling.R|@/data/*|@/inst/*)
            broad_scope=1
            ;;
        @/tests/*)
            broad_scope=1
            ;;
    esac
done < "$scope_file"

LC_ALL=C sort -u "$focus_tmp" > "$focus_file"
rm -f "$focus_tmp" "$changed_file"

{
    printf '\n## Selection notes\n\n'
    if (( broad_scope )); then
        printf -- '- At least one helper, runner, fixture, package-data, installed-resource, or otherwise broad test artifact requires complete-suite coverage.\n'
    fi
    if (( deleted_test_scope )); then
        printf -- '- At least one deleted test remains in candidate scope but was not selected as a runnable focused test; complete-suite coverage is required.\n'
    fi
    if [[ -s "$focus_file" ]]; then
        printf -- '- Focused candidates were selected from changed test files and `tested_by` results. Agent review remains required.\n'
    else
        printf -- '- No focused candidate was selected. An empty focused manifest is allowed only when scope is absent, broad, explicitly full-suite-only, or has no reliable mapping.\n'
    fi
    printf -- '- `tests/spelling.R` and installed-package execution through `tests/testthat.R` are outside this source-level testthat workflow and require separate package-level verification when affected.\n'
} >> "$plan_file"

scope_count="$(wc -l < "$scope_file" | tr -d ' ')"
focus_count="$(wc -l < "$focus_file" | tr -d ' ')"
printf 'SCOPE_COUNT=%s\nFOCUSED_COUNT=%s\n' "$scope_count" "$focus_count"
printf '\n%s\n' "$(sed -n '1,260p' "$plan_file")"
if [[ -s "$focus_file" ]]; then
    printf '\nFOCUSED_TESTS:\n'
    sed -n '1,260p' "$focus_file"
fi
```

## Repair Test-Plan Discovery

```toml donna
id = "repair_test_plan_discovery"
kind = "donna.lib.request_action"
```

Automatic test-plan discovery failed.

Standard output:

```text
{{ donna.lib.task_variable("discovery_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("discovery_stderr") }}
```

Repair only an incorrect discovery rule, a malformed in-scope path, or a Depmesh query discrepancy. Query all configured relations before editing a governed project artifact. Preserve unrelated changes and do not inspect generated manuals.

After repair, `{{ donna.lib.goto("verify_workflow_artifact") }}`. If scope or ownership cannot be determined safely, `{{ donna.lib.goto("record_blocker") }}`.

## Review Test Plan

```toml donna
id = "review_test_plan"
kind = "donna.lib.request_action"
```

Automatic discovery produced a session-local test plan.

Discovery output:

```text
{{ donna.lib.task_variable("discovery_stdout") }}
```

1. Inspect `{{ donna.lib.path("@/.session/donna/run-tests/scope.txt") }}`, `{{ donna.lib.path("@/.session/donna/run-tests/focused-tests.txt") }}`, and `{{ donna.lib.path("@/.session/donna/run-tests/plan.md") }}`. Treat every worktree-derived scope entry as a candidate, not authorization.
2. Query all configured relations separately for every scope artifact and focused test. For a test, inspect its reverse `tests` result; for an R source, inspect its `tested_by` result.
3. Remove every unrelated or pre-existing worktree entry, then add initiating artifacts that are known from the request but absent from the discovered candidates. Use one sorted, unique root-anchored artifact identifier per line.
4. Keep focused entries limited to existing `@/tests/testthat/test-<context>.R` files. A changed test includes itself; changed R sources normally use `tested_by`.
5. Treat empty Depmesh results as unknown, not as proof that no focused test exists. Research consumers of helpers, fixtures, or package data, and record the selection rationale in `plan.md`.
6. It is valid to leave the focused manifest empty when scope is unknown, broad, explicitly full-suite-only, or lacks a reliable focused mapping. The complete suite will still run.
7. Confirm that only entries justified by the initiating request remain. Completing this review makes that pruned scope the maximum authorized repair boundary; dirty worktree state alone MUST NOT authorize a repair.
8. If the plan is accurate, `{{ donna.lib.goto("validate_test_plan") }}`.
9. If the session manifests or rationale need correction before validation, `{{ donna.lib.goto("edit_test_plan") }}`.
10. If safe selection requires unavailable context or a product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Edit Test Plan

```toml donna
id = "edit_test_plan"
kind = "donna.lib.request_action"
```

Edit only the session-local scope, focused-test, and rationale files under `{{ donna.lib.path("@/.session/donna/run-tests") }}`. Query Depmesh separately for every added artifact. Keep entries sorted and unique, do not add generated files or paths outside the project, and do not change package behavior.

After editing, `{{ donna.lib.goto("validate_test_plan") }}`. If the correct scope cannot be established, `{{ donna.lib.goto("record_blocker") }}`.

## Repair Test Plan

```toml donna
id = "repair_test_plan"
kind = "donna.lib.request_action"
```

Test-plan validation output:

```text
{{ donna.lib.task_variable("plan_stdout") }}
```

Test-plan validation diagnostics:

```text
{{ donna.lib.task_variable("plan_stderr") }}
```

1. If a session-local scope, focused-test, or rationale file is incorrect, edit only those files under `{{ donna.lib.path("@/.session/donna/run-tests") }}`. Query Depmesh separately for every added artifact. Keep entries sorted and unique, do not add generated files or paths outside the project, and do not change package behavior. Then `{{ donna.lib.goto("validate_test_plan") }}`.
2. If the validation or discovery rule in this workflow is incorrect, `{{ donna.lib.goto("repair_workflow_artifact") }}`.
3. If the correct scope cannot be established or Depmesh reports an out-of-scope configuration problem, `{{ donna.lib.goto("record_blocker") }}`.

## Validate Test Plan

```toml donna
id = "validate_test_plan"
kind = "donna.lib.run_script"
save_stdout_to = "plan_stdout"
save_stderr_to = "plan_stderr"
goto_on_success = "run_focused_tests"
goto_on_failure = "repair_test_plan"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/run-tests"
scope_file="$plan_dir/scope.txt"
focus_file="$plan_dir/focused-tests.txt"
plan_file="$plan_dir/plan.md"

test -f "$scope_file"
test -f "$focus_file"
test -s "$plan_file"

for manifest in "$scope_file" "$focus_file"; do
    if rg -n '^$|[[:blank:]]+$' "$manifest"; then
        printf 'Manifest contains a blank line or trailing whitespace: %s\n' \
            "$manifest" >&2
        exit 1
    fi
    if ! cmp -s "$manifest" <(LC_ALL=C sort -u "$manifest"); then
        printf 'Manifest must be sorted and unique: %s\n' "$manifest" >&2
        exit 1
    fi
done

while IFS= read -r artifact_id; do
    [[ -n "$artifact_id" ]] || continue
    if [[ ! "$artifact_id" =~ ^@/(R|tests|data|inst)/[^[:space:]]+$ ]] ||
       [[ "$artifact_id" == *"/../"* ]] ||
       [[ "$artifact_id" == "@/man" || "$artifact_id" == "@/man/"* ]]; then
        printf 'Invalid scope artifact identifier: %s\n' "$artifact_id" >&2
        exit 1
    fi

    relative="${artifact_id#@/}"
    unstaged_deletions="$(
        git diff --name-only --diff-filter=D -- "$relative"
    )"
    staged_deletions="$(
        git diff --cached --name-only --diff-filter=D -- "$relative"
    )"
    if [[ ! -e "$relative" ]] &&
       ! printf '%s\n' "$unstaged_deletions" | rg -Fqx -- "$relative" &&
       ! printf '%s\n' "$staged_deletions" | rg -Fqx -- "$relative"; then
        printf 'Scope artifact is neither present nor an unstaged or staged tracked deletion: %s\n' \
            "$artifact_id" >&2
        exit 1
    fi

    query_output="$(depmesh -p llm dependencies "$artifact_id")"
    if printf '%s\n' "$query_output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for scope artifact %s:\n%s\n' \
            "$artifact_id" "$query_output" >&2
        exit 1
    fi
done < "$scope_file"

while IFS= read -r test_id; do
    [[ -n "$test_id" ]] || continue
    if [[ ! "$test_id" =~ ^@/tests/testthat/test-[A-Za-z0-9_-]+\.R$ ]]; then
        printf 'Invalid focused test identifier: %s\n' "$test_id" >&2
        exit 1
    fi
    test_path="${test_id#@/}"
    if [[ ! -f "$test_path" ]]; then
        printf 'Focused test does not exist: %s\n' "$test_id" >&2
        exit 1
    fi
    query_output="$(depmesh -p llm dependencies "$test_id")"
    if printf '%s\n' "$query_output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for focused test %s:\n%s\n' \
            "$test_id" "$query_output" >&2
        exit 1
    fi
done < "$focus_file"

scope_count="$(wc -l < "$scope_file" | tr -d ' ')"
focus_count="$(wc -l < "$focus_file" | tr -d ' ')"
printf 'Validated session test plan with %s scope artifacts and %s focused tests.\n' \
    "$scope_count" "$focus_count"
```

## Run Focused Tests

```toml donna
id = "run_focused_tests"
kind = "donna.lib.run_script"
save_stdout_to = "focused_stdout"
save_stderr_to = "focused_stderr"
goto_on_success = "select_complete_suite_executor"
goto_on_failure = "diagnose_focused_failure"
timeout = 900

[goto_on_code]
"20" = "select_complete_suite_executor"
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

focus_file=".session/donna/run-tests/focused-tests.txt"
test -f "$focus_file"

if [[ ! -s "$focus_file" ]]; then
    printf 'FOCUSED_STATUS=not-run\nFOCUSED_REASON=no reliable focused scope was selected\n'
    exit 20
fi

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected graphics artifact before or after focused tests: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

remove_generated_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -L "$artifact" ]]; then
            printf 'Refusing to remove a graphics-artifact symlink: %s\n' \
                "$artifact" >&2
            return 1
        fi
        if [[ -e "$artifact" ]]; then
            if [[ ! -f "$artifact" ]]; then
                printf 'Refusing to remove a non-regular graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
            if rm -f -- "$artifact"; then
                printf 'REMOVED_TEST_GRAPHICS_ARTIFACT=%s\n' "$artifact"
            else
                printf 'Failed to remove test graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
        fi
    done
    assert_absent_graphics_artifacts
}

cleanup_graphics_artifacts_on_exit() {
    local status=$?
    trap - EXIT
    if ! remove_generated_graphics_artifacts; then
        status=1
    fi
    exit "$status"
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

sha256sum --check .session/donna/run-tests/environment.sha256
assert_absent_graphics_artifacts
trap cleanup_graphics_artifacts_on_exit EXIT
if [[ "$(hash_non_man_worktree)" != \
      "$(cat .session/donna/run-tests/worktree.sha256)" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot before focused-test execution.\n' >&2
    exit 1
fi

pixi run --as-is Rscript --vanilla -e '
options(mc.cores = 1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
if (!inherits(BiocParallel::bpparam(), "SerialParam")) {
    stop("The default BiocParallel backend is not SerialParam")
}
ids <- readLines(".session/donna/run-tests/focused-tests.txt", warn = FALSE)
contexts <- sub("[.]R$", "", sub("^test-", "", basename(ids)))
if (!length(contexts) || any(!grepl("^[A-Za-z0-9_-]+$", contexts))) {
    stop("Focused test contexts are empty or unsafe")
}
filter <- paste0("^(", paste(contexts, collapse = "|"), ")$")
cat("FOCUSED_FILTER=", filter, "\n", sep = "")
results <- testthat::test_local(
    ".",
    filter = filter,
    reporter = "summary",
    load_package = "source",
    stop_on_failure = TRUE,
    stop_on_warning = TRUE
)
summary <- as.data.frame(results)
failed <- sum(summary$failed, na.rm = TRUE)
errors <- sum(summary$error, na.rm = TRUE)
warnings <- sum(summary$warning, na.rm = TRUE)
skipped <- sum(summary$skipped, na.rm = TRUE)
cat(
    "FOCUSED_SUMMARY tests=", nrow(summary),
    " failed=", failed,
    " errors=", errors,
    " warnings=", warnings,
    " skipped=", skipped,
    "\n",
    sep = ""
)
if (!nrow(summary)) {
    stop("The focused filter matched zero tests")
}
if (failed > 0L || errors > 0L || warnings > 0L) {
    stop("Focused testthat results are not successful")
}
if (any(summary$skipped > 0, na.rm = TRUE)) {
    cat("FOCUSED_SKIPPED_TESTS:\n")
    print(summary[summary$skipped > 0, c("file", "test"), drop = FALSE],
          row.names = FALSE)
}
'

sha256sum --check .session/donna/run-tests/environment.sha256
post_focused_failure=0
if ! remove_generated_graphics_artifacts; then
    post_focused_failure=1
fi
if [[ "$(hash_non_man_worktree)" != \
      "$(cat .session/donna/run-tests/worktree.sha256)" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot after focused-test execution.\n' >&2
    post_focused_failure=1
else
    printf 'FOCUSED_GIT_VISIBLE_NON_MAN_WORKTREE=unchanged\n'
fi
if (( post_focused_failure )); then
    exit 1
fi
trap - EXIT
```

## Diagnose Focused Failure

```toml donna
id = "diagnose_focused_failure"
kind = "donna.lib.request_action"
```

Focused test execution failed.

Standard output:

```text
{{ donna.lib.task_variable("focused_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("focused_stderr") }}
```

1. Identify expectation failures, errors, unexpected warnings, package-load failures, timeouts, environment mutations, or graphics artifacts.
2. Query all configured relations separately for every failing test and affected R source.
3. Confirm whether the failure is inside the authorized scope recorded in the session test plan.
4. If the focused selection itself is wrong, `{{ donna.lib.goto("repair_test_plan") }}`.
5. If an in-scope source, test, fixture, or package-data repair is justified, `{{ donna.lib.goto("repair_in_scope") }}`.
6. If the failure exposes an incorrect test-harness or environment-isolation rule in this workflow, `{{ donna.lib.goto("repair_workflow_artifact") }}`.
7. If the failure is pre-existing, out of scope, environment-dependent, or requires generated manuals, installation, network access, or a product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Select Complete-Suite Executor

```toml donna
id = "select_complete_suite_executor"
kind = "donna.lib.request_action"
```

The complete source-level testthat suite remains mandatory.

1. If no explicit resource or authority boundary prevents agent-side execution, `{{ donna.lib.goto("run_complete_suite") }}`.
2. If the maintainer has required the memory-heavy R process to run outside the agent while Donna remains agent-operated, `{{ donna.lib.goto("prepare_manual_complete_suite") }}`.
3. If neither execution path is available, `{{ donna.lib.goto("record_blocker") }}`.

## Prepare Manual Complete Suite

```toml donna
id = "prepare_manual_complete_suite"
kind = "donna.lib.run_script"
save_stdout_to = "manual_prepare_stdout"
save_stderr_to = "manual_prepare_stderr"
goto_on_success = "await_manual_complete_suite"
goto_on_failure = "diagnose_manual_prepare_failure"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/run-tests"

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected graphics artifact before manual-suite preparation: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

validate_environment_manifest() {
    local manifest="$1"
    [[ -f "$manifest" && ! -L "$manifest" ]]
    [[ "$(wc -l < "$manifest")" == "2" ]]
    [[ "$(rg -c '^[0-9a-f]{64}  pixi[.]toml$' "$manifest" || true)" == "1" ]]
    [[ "$(rg -c '^[0-9a-f]{64}  pixi[.]lock$' "$manifest" || true)" == "1" ]]
}

validate_digest_file() {
    local digest_file="$1"
    [[ -f "$digest_file" && ! -L "$digest_file" ]]
    [[ "$(wc -l < "$digest_file")" == "1" ]]
    rg -q '^[0-9a-f]{64}$' "$digest_file"
}

validate_environment_manifest "$plan_dir/environment.sha256"
validate_digest_file "$plan_dir/worktree.sha256"
sha256sum --check "$plan_dir/environment.sha256"
assert_absent_graphics_artifacts
if [[ "$(hash_non_man_worktree)" != "$(cat "$plan_dir/worktree.sha256")" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot before manual-suite preparation.\n' >&2
    exit 1
fi

attempt_dir="$(mktemp -d "$plan_dir/manual-complete.XXXXXXXX")"
attempt_id="$(basename -- "$attempt_dir")"
runner_path="$attempt_dir/run.sh"

cp -- "$plan_dir/environment.sha256" \
    "$attempt_dir/expected-environment.sha256"
cp -- "$plan_dir/worktree.sha256" \
    "$attempt_dir/expected-worktree.sha256"

cat > "$runner_path" <<'RUNNER'
#!/usr/bin/env bash
set -uo pipefail

runner_path="${BASH_SOURCE[0]}"
attempt_dir="$(cd -- "$(dirname -- "$runner_path")" && pwd -P)"
attempt_id="$(basename -- "$attempt_dir")"
project_root="$(git -C "$attempt_dir" rev-parse --show-toplevel)"
log_file="$attempt_dir/combined.log"
evidence_file="$attempt_dir/evidence.tsv"
started_file="$attempt_dir/started"
expected_environment="$attempt_dir/expected-environment.sha256"
expected_worktree="$attempt_dir/expected-worktree.sha256"

cd -- "$project_root"

if ! (
    set -o noclobber
    printf '%s\n' "$attempt_id" > "$started_file"
) 2>/dev/null; then
    printf 'This manual-suite attempt has already started: %s\n' \
        "$attempt_id" >&2
    exit 2
fi

for evidence_input in \
    "$attempt_dir/runner.sha256" \
    "$expected_environment" \
    "$expected_worktree"
do
    if [[ ! -f "$evidence_input" || -L "$evidence_input" ]]; then
        printf 'Invalid manual-suite evidence input: %s\n' \
            "$evidence_input" >&2
        exit 2
    fi
done
runner_expected="$(cat "$attempt_dir/runner.sha256")"
runner_actual="$(sha256sum -- "$runner_path" | awk '{ print $1 }')"
if [[ ! "$runner_expected" =~ ^[0-9a-f]{64}$ ||
      "$(wc -l < "$attempt_dir/runner.sha256")" != "1" ||
      "$runner_actual" != "$runner_expected" ]]; then
    printf 'Manual-suite runner checksum verification failed.\n' >&2
    exit 2
fi
if [[ "$(wc -l < "$expected_environment")" != "2" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]toml$' \
          "$expected_environment" || true)" != "1" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]lock$' \
          "$expected_environment" || true)" != "1" ]]; then
    printf 'Invalid expected-environment manifest.\n' >&2
    exit 2
fi
if [[ "$(wc -l < "$expected_worktree")" != "1" ]] ||
   ! rg -q '^[0-9a-f]{64}$' "$expected_worktree"; then
    printf 'Invalid expected-worktree digest.\n' >&2
    exit 2
fi

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected graphics artifact before or after complete tests: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

remove_generated_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -L "$artifact" ]]; then
            printf 'Refusing to remove a graphics-artifact symlink: %s\n' \
                "$artifact" >&2
            return 1
        fi
        if [[ -e "$artifact" ]]; then
            if [[ ! -f "$artifact" ]]; then
                printf 'Refusing to remove a non-regular graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
            if rm -f -- "$artifact"; then
                printf 'REMOVED_TEST_GRAPHICS_ARTIFACT=%s\n' "$artifact"
            else
                printf 'Failed to remove test graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
        fi
    done
    assert_absent_graphics_artifacts
}

cleanup_graphics_artifacts_on_exit() {
    local status=$?
    trap - EXIT
    if ! remove_generated_graphics_artifacts; then
        status=1
    fi
    exit "$status"
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

run_complete_suite() {
    set -euo pipefail

    sha256sum --check "$expected_environment"
    assert_absent_graphics_artifacts
    trap cleanup_graphics_artifacts_on_exit EXIT
    if [[ "$(hash_non_man_worktree)" != "$(cat "$expected_worktree")" ]]; then
        printf 'The Git-visible non-man worktree differs from its pre-test snapshot before maintainer-executed complete-suite execution.\n' >&2
        exit 1
    fi

    printf 'COMPLETE_EXECUTOR=maintainer\n'
    printf 'MANUAL_COMPLETE_ATTEMPT=%s\n' "$attempt_id"
    pixi run --as-is Rscript --vanilla -e '
options(mc.cores = 1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
if (!inherits(BiocParallel::bpparam(), "SerialParam")) {
    stop("The default BiocParallel backend is not SerialParam")
}
results <- testthat::test_local(
    ".",
    reporter = "summary",
    load_package = "source",
    stop_on_failure = TRUE,
    stop_on_warning = TRUE
)
summary <- as.data.frame(results)
failed <- sum(summary$failed, na.rm = TRUE)
errors <- sum(summary$error, na.rm = TRUE)
warnings <- sum(summary$warning, na.rm = TRUE)
skipped <- sum(summary$skipped, na.rm = TRUE)
cat(
    "COMPLETE_SUMMARY tests=", nrow(summary),
    " failed=", failed,
    " errors=", errors,
    " warnings=", warnings,
    " skipped=", skipped,
    "\n",
    sep = ""
)
if (!nrow(summary)) {
    stop("The complete suite executed zero tests")
}
if (failed > 0L || errors > 0L || warnings > 0L) {
    stop("Complete testthat results are not successful")
}
if (any(summary$skipped > 0, na.rm = TRUE)) {
    cat("COMPLETE_SKIPPED_TESTS:\n")
    print(summary[summary$skipped > 0, c("file", "test"), drop = FALSE],
          row.names = FALSE)
}
'

    sha256sum --check "$expected_environment"
    post_complete_failure=0
    if ! remove_generated_graphics_artifacts; then
        post_complete_failure=1
    fi
    if [[ "$(hash_non_man_worktree)" != "$(cat "$expected_worktree")" ]]; then
        printf 'The Git-visible non-man worktree differs from its pre-test snapshot after maintainer-executed complete-suite execution.\n' >&2
        post_complete_failure=1
    else
        printf 'COMPLETE_GIT_VISIBLE_NON_MAN_WORKTREE=unchanged\n'
    fi
    if (( post_complete_failure )); then
        exit 1
    fi
    trap - EXIT
}

: > "$log_file"
set +e
(
    run_complete_suite
) 2>&1 | tee "$log_file"
pipeline_status=("${PIPESTATUS[@]}")
suite_status="${pipeline_status[0]}"
log_status="${pipeline_status[1]}"
set -e

environment_status=0
if ! sha256sum --check "$expected_environment" >/dev/null 2>&1; then
    environment_status=1
fi
worktree_status=0
if [[ "$(hash_non_man_worktree)" != "$(cat "$expected_worktree")" ]]; then
    worktree_status=1
fi
graphics_status=0
for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
    if [[ -e "$artifact" || -L "$artifact" ]]; then
        graphics_status=1
    fi
done

overall_status="$suite_status"
if (( log_status || environment_status || worktree_status || graphics_status )); then
    overall_status=1
fi

log_sha256="$(sha256sum -- "$log_file" | awk '{ print $1 }')"
evidence_tmp="$attempt_dir/evidence.tsv.tmp"
{
    printf 'schema=1\n'
    printf 'attempt=%s\n' "$attempt_id"
    printf 'suite_status=%s\n' "$suite_status"
    printf 'log_status=%s\n' "$log_status"
    printf 'overall_status=%s\n' "$overall_status"
    printf 'environment_status=%s\n' "$environment_status"
    printf 'worktree_status=%s\n' "$worktree_status"
    printf 'graphics_status=%s\n' "$graphics_status"
    printf 'log_sha256=%s\n' "$log_sha256"
    printf 'complete=1\n'
} > "$evidence_tmp"
mv -- "$evidence_tmp" "$evidence_file"

printf 'MANUAL_COMPLETE_ATTEMPT=%s\n' "$attempt_id"
printf 'MANUAL_COMPLETE_STATUS=%s\n' "$overall_status"
exit "$overall_status"
RUNNER

chmod 700 "$runner_path"
runner_sha256="$(sha256sum -- "$runner_path" | awk '{ print $1 }')"
environment_manifest_sha256="$(
    sha256sum -- "$attempt_dir/expected-environment.sha256" |
        awk '{ print $1 }'
)"
worktree_digest_sha256="$(
    sha256sum -- "$attempt_dir/expected-worktree.sha256" |
        awk '{ print $1 }'
)"
printf '%s\n' "$runner_sha256" > "$attempt_dir/runner.sha256"

printf 'MANUAL_COMPLETE_ATTEMPT=%s\n' "$attempt_id"
printf 'MANUAL_RUNNER_SHA256=%s\n' "$runner_sha256"
printf 'MANUAL_ENVIRONMENT_MANIFEST_SHA256=%s\n' \
    "$environment_manifest_sha256"
printf 'MANUAL_WORKTREE_DIGEST_SHA256=%s\n' "$worktree_digest_sha256"
printf 'MANUAL_COMPLETE_COMMAND=bash %q\n' "$runner_path"
```

## Diagnose Manual Preparation Failure

```toml donna
id = "diagnose_manual_prepare_failure"
kind = "donna.lib.request_action"
```

Manual complete-suite preparation failed.

Standard output:

```text
{{ donna.lib.task_variable("manual_prepare_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("manual_prepare_stderr") }}
```

1. If the workflow or its session-evidence rules are incorrect, `{{ donna.lib.goto("repair_workflow_artifact") }}`.
2. If the worktree or Pixi snapshot changed intentionally, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
3. If preparation is blocked by a pre-existing graphics artifact, unsafe filesystem object, or unavailable local tool, `{{ donna.lib.goto("record_blocker") }}`.

## Await Manual Complete Suite

```toml donna
id = "await_manual_complete_suite"
kind = "donna.lib.request_action"
```

Donna prepared a unique maintainer-executed complete-suite attempt:

```text
{{ donna.lib.task_variable("manual_prepare_stdout") }}
```

Preparation diagnostics:

```text
{{ donna.lib.task_variable("manual_prepare_stderr") }}
```

1. Copy the value after `MANUAL_COMPLETE_COMMAND=` exactly and give it to the maintainer in a shell code block. The maintainer runs only that command from the project root; the maintainer does not run Donna.
2. The agent MUST NOT execute the prepared command. Leave this action request pending while the maintainer runs it and returns the terminal output. Do not run another test process or create either declared `Rplots.pdf` path concurrently with this attempt.
3. Confirm the returned `MANUAL_COMPLETE_ATTEMPT` matches the prepared attempt. Use the session-local evidence as the verification source; do not reconstruct evidence from pasted output.
4. After the maintainer reports that the command finished, `{{ donna.lib.goto("verify_manual_complete_suite") }}`.
5. If the worktree changed before the command could run, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
6. If the maintainer cannot run the command, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Manual Complete Suite

```toml donna
id = "verify_manual_complete_suite"
kind = "donna.lib.run_script"
save_stdout_to = "complete_stdout"
save_stderr_to = "complete_stderr"
goto_on_success = "verify_runtime_state"
goto_on_failure = "diagnose_complete_failure"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/run-tests"
prepared_output="$(cat <<'DONNA_MANUAL_PREPARE_OUTPUT'
{{ donna.lib.task_variable("manual_prepare_stdout") }}
DONNA_MANUAL_PREPARE_OUTPUT
)"

prepare_value() {
    local key="$1"
    local count value
    count="$(
        printf '%s\n' "$prepared_output" |
            awk -F= -v key="$key" \
                '$1 == key { count++ } END { print count + 0 }'
    )"
    if [[ "$count" != "1" ]]; then
        printf 'Manual preparation key %s occurs %s times.\n' \
            "$key" "$count" >&2
        return 1
    fi
    value="$(
        printf '%s\n' "$prepared_output" |
            awk -F= -v key="$key" \
                '$1 == key { print substr($0, length(key) + 2) }'
    )"
    printf '%s\n' "$value"
}

attempt_id="$(prepare_value MANUAL_COMPLETE_ATTEMPT)"
prepared_runner_sha256="$(prepare_value MANUAL_RUNNER_SHA256)"
prepared_environment_manifest_sha256="$(
    prepare_value MANUAL_ENVIRONMENT_MANIFEST_SHA256
)"
prepared_worktree_digest_sha256="$(
    prepare_value MANUAL_WORKTREE_DIGEST_SHA256
)"
if [[ ! "$attempt_id" =~ ^manual-complete\.[A-Za-z0-9]{8}$ ]]; then
    printf 'Invalid manual complete-suite attempt identifier: %s\n' \
        "$attempt_id" >&2
    exit 1
fi
for prepared_digest in \
    "$prepared_runner_sha256" \
    "$prepared_environment_manifest_sha256" \
    "$prepared_worktree_digest_sha256"
do
    if [[ ! "$prepared_digest" =~ ^[0-9a-f]{64}$ ]]; then
        printf 'Invalid manual preparation digest: %s\n' \
            "$prepared_digest" >&2
        exit 1
    fi
done

attempt_dir="$plan_dir/$attempt_id"
if [[ ! -d "$attempt_dir" || -L "$attempt_dir" ]]; then
    printf 'Invalid manual complete-suite attempt directory: %s\n' \
        "$attempt_dir" >&2
    exit 1
fi

for required_file in \
    run.sh \
    runner.sha256 \
    expected-environment.sha256 \
    expected-worktree.sha256 \
    started \
    combined.log \
    evidence.tsv
do
    required_path="$attempt_dir/$required_file"
    if [[ ! -f "$required_path" || -L "$required_path" ]]; then
        printf 'Missing or unsafe manual-suite evidence file: %s\n' \
            "$required_path" >&2
        exit 1
    fi
done

if [[ "$(wc -l < "$attempt_dir/started")" != "1" ||
      "$(cat "$attempt_dir/started")" != "$attempt_id" ]]; then
    printf 'Manual-suite start marker is invalid.\n' >&2
    exit 1
fi

runner_expected="$(cat "$attempt_dir/runner.sha256")"
runner_actual="$(sha256sum -- "$attempt_dir/run.sh" | awk '{ print $1 }')"
if [[ ! "$runner_expected" =~ ^[0-9a-f]{64}$ ||
      "$(wc -l < "$attempt_dir/runner.sha256")" != "1" ||
      "$runner_expected" != "$prepared_runner_sha256" ||
      "$runner_actual" != "$runner_expected" ]]; then
    printf 'Manual-suite runner checksum verification failed.\n' >&2
    exit 1
fi
if [[ "$(wc -l < "$attempt_dir/expected-environment.sha256")" != "2" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]toml$' \
          "$attempt_dir/expected-environment.sha256" || true)" != "1" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]lock$' \
          "$attempt_dir/expected-environment.sha256" || true)" != "1" ]]; then
    printf 'Invalid expected-environment manifest.\n' >&2
    exit 1
fi
if [[ "$(wc -l < "$attempt_dir/expected-worktree.sha256")" != "1" ]] ||
   ! rg -q '^[0-9a-f]{64}$' \
       "$attempt_dir/expected-worktree.sha256"; then
    printf 'Invalid expected-worktree digest.\n' >&2
    exit 1
fi
if [[ "$(sha256sum -- "$attempt_dir/expected-environment.sha256" |
          awk '{ print $1 }')" != "$prepared_environment_manifest_sha256" ||
      "$(sha256sum -- "$attempt_dir/expected-worktree.sha256" |
          awk '{ print $1 }')" != "$prepared_worktree_digest_sha256" ]]; then
    printf 'Manual-suite snapshot seal does not match Donna preparation.\n' >&2
    exit 1
fi
if [[ ! -f "$plan_dir/environment.sha256" ||
      -L "$plan_dir/environment.sha256" ||
      "$(wc -l < "$plan_dir/environment.sha256")" != "2" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]toml$' \
          "$plan_dir/environment.sha256" || true)" != "1" ||
      "$(rg -c '^[0-9a-f]{64}  pixi[.]lock$' \
          "$plan_dir/environment.sha256" || true)" != "1" ||
      ! -f "$plan_dir/worktree.sha256" ||
      -L "$plan_dir/worktree.sha256" ||
      "$(wc -l < "$plan_dir/worktree.sha256")" != "1" ]] ||
   ! rg -q '^[0-9a-f]{64}$' "$plan_dir/worktree.sha256"; then
    printf 'Current Donna snapshot files are missing or unsafe.\n' >&2
    exit 1
fi
if ! cmp -s -- "$attempt_dir/expected-environment.sha256" \
        "$plan_dir/environment.sha256" ||
   ! cmp -s -- "$attempt_dir/expected-worktree.sha256" \
        "$plan_dir/worktree.sha256"; then
    printf 'Manual-suite snapshots differ from the current Donna plan.\n' >&2
    exit 1
fi

cat "$attempt_dir/combined.log"

evidence_file="$attempt_dir/evidence.tsv"
if (( $(wc -l < "$evidence_file") != 10 )); then
    printf 'Manual-suite evidence must contain exactly ten records.\n' >&2
    exit 1
fi

evidence_value() {
    local key="$1"
    local count value
    count="$(awk -F= -v key="$key" '$1 == key { count++ } END { print count + 0 }' \
        "$evidence_file")"
    if [[ "$count" != "1" ]]; then
        printf 'Manual-suite evidence key %s occurs %s times.\n' \
            "$key" "$count" >&2
        return 1
    fi
    value="$(awk -F= -v key="$key" \
        '$1 == key { print substr($0, length(key) + 2) }' "$evidence_file")"
    printf '%s\n' "$value"
}

schema="$(evidence_value schema)"
evidence_attempt="$(evidence_value attempt)"
suite_status="$(evidence_value suite_status)"
log_status="$(evidence_value log_status)"
overall_status="$(evidence_value overall_status)"
environment_status="$(evidence_value environment_status)"
worktree_status="$(evidence_value worktree_status)"
graphics_status="$(evidence_value graphics_status)"
log_sha256="$(evidence_value log_sha256)"
complete="$(evidence_value complete)"

[[ "$schema" == "1" ]]
[[ "$evidence_attempt" == "$attempt_id" ]]
[[ "$suite_status" == "0" ]]
[[ "$log_status" == "0" ]]
[[ "$overall_status" == "0" ]]
[[ "$environment_status" == "0" ]]
[[ "$worktree_status" == "0" ]]
[[ "$graphics_status" == "0" ]]
[[ "$complete" == "1" ]]
[[ "$log_sha256" =~ ^[0-9a-f]{64}$ ]]
[[ "$log_sha256" == \
   "$(sha256sum -- "$attempt_dir/combined.log" | awk '{ print $1 }')" ]]

summary_count="$(
    rg -c \
        '^COMPLETE_SUMMARY tests=[1-9][0-9]* failed=[0-9]+ errors=[0-9]+ warnings=[0-9]+ skipped=[0-9]+$' \
        "$attempt_dir/combined.log" ||
        true
)"
if [[ "$summary_count" != "1" ]]; then
    printf 'Expected exactly one machine-readable complete-suite summary; found %s.\n' \
        "$summary_count" >&2
    exit 1
fi
if [[ "$(rg -Fxc -- 'COMPLETE_EXECUTOR=maintainer' \
          "$attempt_dir/combined.log" || true)" != "1" ||
      "$(rg -Fxc -- "MANUAL_COMPLETE_ATTEMPT=$attempt_id" \
          "$attempt_dir/combined.log" || true)" != "1" ]]; then
    printf 'Manual-suite executor or attempt marker is missing or duplicated.\n' >&2
    exit 1
fi
if ! rg -q \
    '^COMPLETE_SUMMARY tests=[1-9][0-9]* failed=0 errors=0 warnings=0 skipped=[0-9]+$' \
    "$attempt_dir/combined.log"
then
    printf 'Complete-suite summary reports failures, errors, or warnings.\n' >&2
    exit 1
fi

sha256sum --check "$attempt_dir/expected-environment.sha256"

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

if [[ "$(hash_non_man_worktree)" != \
      "$(cat "$attempt_dir/expected-worktree.sha256")" ]]; then
    printf 'The Git-visible non-man worktree differs from the manual-suite snapshot during evidence verification.\n' >&2
    exit 1
fi
for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
    if [[ -e "$artifact" || -L "$artifact" ]]; then
        printf 'Unexpected post-test graphics artifact: %s\n' \
            "$artifact" >&2
        exit 1
    fi
done

printf 'MANUAL_COMPLETE_EVIDENCE=verified\n'
printf 'MANUAL_COMPLETE_ATTEMPT=%s\n' "$attempt_id"
```

## Run Complete Suite

```toml donna
id = "run_complete_suite"
kind = "donna.lib.run_script"
save_stdout_to = "complete_stdout"
save_stderr_to = "complete_stderr"
goto_on_success = "verify_runtime_state"
goto_on_failure = "diagnose_complete_failure"
timeout = 1800
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected graphics artifact before or after complete tests: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

remove_generated_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -L "$artifact" ]]; then
            printf 'Refusing to remove a graphics-artifact symlink: %s\n' \
                "$artifact" >&2
            return 1
        fi
        if [[ -e "$artifact" ]]; then
            if [[ ! -f "$artifact" ]]; then
                printf 'Refusing to remove a non-regular graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
            if rm -f -- "$artifact"; then
                printf 'REMOVED_TEST_GRAPHICS_ARTIFACT=%s\n' "$artifact"
            else
                printf 'Failed to remove test graphics artifact: %s\n' \
                    "$artifact" >&2
                return 1
            fi
        fi
    done
    assert_absent_graphics_artifacts
}

cleanup_graphics_artifacts_on_exit() {
    local status=$?
    trap - EXIT
    if ! remove_generated_graphics_artifacts; then
        status=1
    fi
    exit "$status"
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

sha256sum --check .session/donna/run-tests/environment.sha256
assert_absent_graphics_artifacts
trap cleanup_graphics_artifacts_on_exit EXIT
if [[ "$(hash_non_man_worktree)" != \
      "$(cat .session/donna/run-tests/worktree.sha256)" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot before complete-suite execution.\n' >&2
    exit 1
fi

printf 'COMPLETE_EXECUTOR=donna\n'
pixi run --as-is Rscript --vanilla -e '
options(mc.cores = 1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
if (!inherits(BiocParallel::bpparam(), "SerialParam")) {
    stop("The default BiocParallel backend is not SerialParam")
}
results <- testthat::test_local(
    ".",
    reporter = "summary",
    load_package = "source",
    stop_on_failure = TRUE,
    stop_on_warning = TRUE
)
summary <- as.data.frame(results)
failed <- sum(summary$failed, na.rm = TRUE)
errors <- sum(summary$error, na.rm = TRUE)
warnings <- sum(summary$warning, na.rm = TRUE)
skipped <- sum(summary$skipped, na.rm = TRUE)
cat(
    "COMPLETE_SUMMARY tests=", nrow(summary),
    " failed=", failed,
    " errors=", errors,
    " warnings=", warnings,
    " skipped=", skipped,
    "\n",
    sep = ""
)
if (!nrow(summary)) {
    stop("The complete suite executed zero tests")
}
if (failed > 0L || errors > 0L || warnings > 0L) {
    stop("Complete testthat results are not successful")
}
if (any(summary$skipped > 0, na.rm = TRUE)) {
    cat("COMPLETE_SKIPPED_TESTS:\n")
    print(summary[summary$skipped > 0, c("file", "test"), drop = FALSE],
          row.names = FALSE)
}
'

sha256sum --check .session/donna/run-tests/environment.sha256
post_complete_failure=0
if ! remove_generated_graphics_artifacts; then
    post_complete_failure=1
fi
if [[ "$(hash_non_man_worktree)" != \
      "$(cat .session/donna/run-tests/worktree.sha256)" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot after complete-suite execution.\n' >&2
    post_complete_failure=1
else
    printf 'COMPLETE_GIT_VISIBLE_NON_MAN_WORKTREE=unchanged\n'
fi
if (( post_complete_failure )); then
    exit 1
fi
trap - EXIT
```

## Diagnose Complete-Suite Failure

```toml donna
id = "diagnose_complete_failure"
kind = "donna.lib.request_action"
```

Complete testthat execution failed.

Standard output:

```text
{{ donna.lib.task_variable("complete_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("complete_stderr") }}
```

1. Classify expectation failures, errors, unexpected warnings, package-load failures, timeouts, resource limitations, environment mutations, or graphics artifacts.
2. Query all configured relations separately for every failing test and affected R source.
3. Compare the failure with the authorized scope and focused result. Do not treat the testing request as permission to repair unrelated package behavior.
4. If a maintainer-executed attempt is missing, interrupted, or malformed without changing the governed worktree, `{{ donna.lib.goto("prepare_manual_complete_suite") }}`.
5. If the saved environment or worktree hash no longer matches because an intentional change occurred, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
6. If an in-scope source, test, fixture, or package-data repair is justified, `{{ donna.lib.goto("repair_in_scope") }}`.
7. If the failure exposes an incorrect test-harness, manual-evidence, or environment-isolation rule in this workflow, `{{ donna.lib.goto("repair_workflow_artifact") }}`.
8. If the failure is pre-existing, out of scope, environment-dependent, or requires generated manuals, installation, network access, or a product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Repair In Scope

```toml donna
id = "repair_in_scope"
kind = "donna.lib.request_action"
```

1. Read the session scope and plan, then query all configured relations separately for every artifact to be edited.
2. Repair only the authorized R source, test, fixture, or package-data behavior demonstrated by the failure.
3. Preserve public API documentation and registration conventions in their editable R sources, and add or update focused testthat coverage when R behavior changes.
4. Do not edit or inspect generated manuals, install dependencies, use the network, run package checks, stage files, create commits, or broaden the product contract.
5. After repair, reread every changed artifact and `{{ donna.lib.goto("verify_workflow_artifact") }}` so scope discovery, focused tests, and the complete suite all repeat.
6. If a correct repair exceeds the authorized scope or needs developer direction, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Runtime State

```toml donna
id = "verify_runtime_state"
kind = "donna.lib.run_script"
save_stdout_to = "runtime_stdout"
save_stderr_to = "runtime_stderr"
goto_on_success = "review_test_results"
goto_on_failure = "repair_runtime_state"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

assert_absent_graphics_artifacts() {
    local artifact
    for artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
        if [[ -e "$artifact" || -L "$artifact" ]]; then
            printf 'Unexpected post-test graphics artifact: %s\n' \
                "$artifact" >&2
            return 1
        fi
    done
}

assert_graphics_artifacts_not_ignored() {
    local ignore_output ignore_status
    set +e
    ignore_output="$(
        git check-ignore --no-index -v -- \
            Rplots.pdf tests/testthat/Rplots.pdf 2>&1
    )"
    ignore_status=$?
    set -e
    case "$ignore_status" in
        0)
            printf 'Test graphics output is hidden by an ignore rule:\n%s\n' \
                "$ignore_output" >&2
            return 1
            ;;
        1)
            return 0
            ;;
        *)
            printf 'Unable to verify test-graphics ignore state:\n%s\n' \
                "$ignore_output" >&2
            return "$ignore_status"
            ;;
    esac
}

hash_non_man_worktree() {
    local digest link_target untracked_file
    if ! untracked_file="$(
        mktemp /tmp/probatch-run-tests-untracked.XXXXXX
    )"; then
        return 1
    fi
    if ! git ls-files --others --exclude-standard -z -- \
        . ':(top,exclude)man/**' > "$untracked_file"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! digest="$(
        {
            printf 'UNSTAGED\0'
            git diff --no-ext-diff --no-textconv --binary --full-index -- \
                . ':(top,exclude)man/**' ||
                exit 1
            printf 'STAGED\0'
            git diff --cached --no-ext-diff --no-textconv --binary \
                --full-index -- . ':(top,exclude)man/**' ||
                exit 1
            while IFS= read -r -d '' path; do
                case "$path" in
                    Rplots.pdf|tests/testthat/Rplots.pdf) continue ;;
                esac
                printf 'UNTRACKED\0%s\0' "$path"
                if [[ -L "$path" ]]; then
                    if ! link_target="$(readlink -- "$path")"; then
                        exit 1
                    fi
                    printf 'SYMLINK\0%s\0' "$link_target"
                elif [[ -f "$path" ]]; then
                    printf 'FILE\0'
                    sha256sum -- "$path" || exit 1
                else
                    printf 'OTHER\0'
                fi
            done < "$untracked_file"
        } |
            sha256sum |
            awk '{ print $1 }'
    )"; then
        rm -f -- "$untracked_file"
        return 1
    fi
    if ! rm -f -- "$untracked_file"; then
        return 1
    fi
    printf '%s\n' "$digest"
}

sha256sum --check .session/donna/run-tests/environment.sha256
assert_absent_graphics_artifacts
assert_graphics_artifacts_not_ignored
if [[ "$(hash_non_man_worktree)" != \
      "$(cat .session/donna/run-tests/worktree.sha256)" ]]; then
    printf 'The Git-visible non-man worktree differs from its pre-test snapshot during runtime verification.\n' >&2
    exit 1
fi

r_output_artifacts="$(
    find tests -type f \
        \( -name '*.Rout' -o -name '*.Rout.fail' -o -name '*.Rout.save' \) \
        -print
)"
if [[ -n "$r_output_artifacts" ]]; then
    printf 'Test execution left R output artifacts under tests/:\n%s\n' \
        "$r_output_artifacts" >&2
    exit 1
fi

git diff --check -- . ':(top,exclude)man/**'
git diff --cached --check -- . ':(top,exclude)man/**'

depmesh -p llm dependencies @/workflows/run-tests.donna.md

printf 'NON_MAN_WORKTREE_STATUS:\n'
git status --short --untracked-files=all -- \
    . ':(top,exclude)man/**'

printf 'Test execution preserved the Git-visible non-man worktree state and Pixi metadata and left no known graphics or R output artifacts.\n'
```

## Repair Runtime State

```toml donna
id = "repair_runtime_state"
kind = "donna.lib.request_action"
```

Post-test runtime verification failed.

Standard output:

```text
{{ donna.lib.task_variable("runtime_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("runtime_stderr") }}
```

1. If an in-scope package or test defect produced the artifact, `{{ donna.lib.goto("repair_in_scope") }}`.
2. If this workflow, its catalog/guidance synchronization, or its verification expectation is wrong, `{{ donna.lib.goto("repair_workflow_artifact") }}`.
3. Do not hide generated artifacts with ignore rules or repair unrelated user changes.
4. If cleanup or repair requires broader authority, generated-manual inspection, installation, network access, or a product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Review Test Results

```toml donna
id = "review_test_results"
kind = "donna.lib.request_action"
```

Environment verification:

```text
{{ donna.lib.task_variable("environment_stdout") }}
```

Environment diagnostics:

```text
{{ donna.lib.task_variable("environment_stderr") }}
```

Focused execution output:

```text
{{ donna.lib.task_variable("focused_stdout") }}
```

Focused execution diagnostics:

```text
{{ donna.lib.task_variable("focused_stderr") }}
```

Complete-suite output:

```text
{{ donna.lib.task_variable("complete_stdout") }}
```

Complete-suite diagnostics:

```text
{{ donna.lib.task_variable("complete_stderr") }}
```

Runtime verification:

```text
{{ donna.lib.task_variable("runtime_stdout") }}
```

1. Confirm the existing Pixi environment was used without installation and with `pixi.toml` and `pixi.lock` unchanged. Treat BiocParallel as an existing test-environment prerequisite, not a dependency this workflow may provision.
2. Confirm focused tests either passed in their own clean R process or were explicitly not run because no reliable focused scope was selected.
3. Confirm the complete testthat suite passed in a separate clean R process with zero failures, errors, and unexpected test warnings. If the maintainer executed Donna's prepared command, confirm `MANUAL_COMPLETE_EVIDENCE=verified` and the matching attempt identifier; otherwise confirm `COMPLETE_EXECUTOR=donna`.
4. Review every reported skip and its reason. A skip that removes the behavior under focused verification is an environment limitation or blocker, not a focused pass.
5. Distinguish benign optional or locale skips from missing declared Imports or other broken required coverage.
6. Report any process-level diagnostics explicitly, including restricted-port or timezone diagnostics that did not become test failures. Confirm that only newly generated regular `Rplots.pdf` files at the two declared test-output locations were removed. Do not claim coverage for `tests/spelling.R`, installed-package execution through `tests/testthat.R`, vignettes, package checks, or generated documentation.
7. Inspect the exact in-scope diff and untracked workflow without modifying the Git index. Preserve unrelated changes and confirm no generated or cache artifact entered scope.
8. While Donna is awaiting this action and no workflow operation is executing, run `donna -p llm validate @/workflows/run-tests.donna.md`, render it in view and analysis modes, and run `donna -p llm list`.
9. Confirm the catalog, agent guidance, scalable Depmesh rules, test plan, clean-process commands, serial backend registration, and repair loops form one coherent implementation.
10. If all test and workflow evidence is sufficient, `{{ donna.lib.goto("finish") }}`.
11. If an in-scope repair is needed, `{{ donna.lib.goto("repair_in_scope") }}`.
12. If the complete suite must be rerun because evidence is incomplete, `{{ donna.lib.goto("select_complete_suite_executor") }}`.
13. If a skip, diagnostic, or failure prevents the promised verification, `{{ donna.lib.goto("record_blocker") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Focused tests selected from the known scope passed or were explicitly unnecessary, and the complete source-level testthat suite passed in a separate clean R process using a deterministic serial BiocParallel backend. The complete process may be Donna-executed or maintainer-executed through Donna's unique prepared command and verified session evidence. Report the selected scope and focused tests, executor and manual-attempt identifier when applicable, complete-suite summary, every skip and environment diagnostic, removal of newly generated regular `Rplots.pdf` test output, Git-visible non-man worktree preservation, runtime-artifact checks, and residual verification intentionally left to other workflows.

If this was the first representative execution and the catalog still says `implementation in progress`, mark it `implemented` only after this finish. Then rerun this workflow in the final catalog state and repeat specific/all workflow validation, view and analysis rendering, workflow discovery, bidirectional governance queries, non-man diff checks, and the complete testthat suite before reporting completion.

## Record Blocker

```toml donna
id = "record_blocker"
kind = "donna.lib.request_action"
```

The workflow has reached a failure, limitation, conflict, or authority boundary that it cannot resolve directly.

1. Preserve the exact failing output, selected scope, completed gates, and decision or external action required.
2. Determine whether a remaining cataloged workflow owns the resolution. If it does, name that workflow and the relevant stage; do not duplicate its work or broaden the current scope.
3. If no remaining workflow resolves the finding, read the workflow-report requirements, query all configured relations for every report artifact to be edited, and create or update an open-questions Markdown report under `{{ donna.lib.path("@/workflows/reports") }}` before continuing. Reports MUST NOT use the `.donna.md` extension or enter the workflow catalog.
4. Preserve unrelated changes and the Git index. Do not inspect generated manuals, change package behavior, install dependencies, use the network, run package-wide checks, stage files, create commits, or reset the Donna session.
5. After recording or routing the blocker, `{{ donna.lib.goto("blocked") }}`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The focused-package-test workflow cannot complete its promised verification. Report completed gates, selected scope, focused results, the exact failure or skip, whether it is in scope, any environment or authority limitation, and the decision or external action required. Do not inspect generated manuals, install packages, use the network, run package-wide checks, stage files, create commits, or reset the Donna session.
