# Verify File Relations

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Verify that Depmesh implements the complete bidirectional, project-local file-relation contract for specifications, R sources, tests, workflows, reports, configuration, and helpers.

This workflow MUST NOT access `man/`, use the network, change package behavior, stage files, create commits, or reset the Donna session.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/meta/general.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}`.
2. Read the built-in Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations for `{{ donna.lib.path("@/depmesh.toml") }}`, `{{ donna.lib.path("@/donna.toml") }}`, `{{ donna.lib.path("@/AGENTS.md") }}`, this workflow, and every existing helper under `{{ donna.lib.path("@/bin/depemesh") }}`.
3. Inspect the current worktree without modifying the Git index. Preserve unrelated user changes.
4. Confirm that no unrelated Donna task or action request is pending. Do not reset the session.
5. Confirm that this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna execution process was active.
6. Do not read, enumerate, compare, or otherwise inspect any path under `man/`.
7. If the relation contract and current scope are understood, `{{ donna.lib.goto("verify_relation_definitions") }}`.
8. If the specifications conflict or required context is unavailable, `{{ donna.lib.goto("blocked") }}`.

## Verify Relation Definitions

```toml donna
id = "verify_relation_definitions"
kind = "donna.lib.run_script"
save_stdout_to = "definitions_stdout"
save_stderr_to = "definitions_stderr"
goto_on_success = "verify_source_test_mappings"
goto_on_failure = "repair_relation_definitions"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

relation_output="$(depmesh -p llm relations)"
actual_ids="$(
    printf '%s\n' "$relation_output" |
        sed -n 's/^## \([a-z_][a-z_]*\)$/\1/p'
)"
expected_ids=$'governed_by\ngoverns\ntested_by\ntests'

if [[ "$actual_ids" != "$expected_ids" ]]; then
    printf 'Unexpected relation definitions.\nExpected:\n%s\nActual:\n%s\n' \
        "$expected_ids" "$actual_ids" >&2
    exit 1
fi

required_descriptions=(
    "Specifications whose requirements apply to the artifact."
    "Project artifacts governed by the specification."
    "Testthat files that verify the R source file."
    "R source files verified by the testthat file."
)

for description in "${required_descriptions[@]}"; do
    if [[ "$relation_output" != *"$description"* ]]; then
        printf 'Missing or incorrect relation description: %s\n' "$description" >&2
        exit 1
    fi
done

printf '%s\n' "$relation_output"
```

## Repair Relation Definitions

```toml donna
id = "repair_relation_definitions"
kind = "donna.lib.request_action"
```

The relation-definition gate failed.

Standard output:

```text
{{ donna.lib.task_variable("definitions_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("definitions_stderr") }}
```

Query all relations for every artifact to be edited, then repair only the relation definitions, their governing specification, or this verifier. Preserve unrelated changes and do not inspect `man/`.

After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`. If repair requires a product decision or conflicting specification change, `{{ donna.lib.goto("blocked") }}`.

## Verify Source and Test Mappings

```toml donna
id = "verify_source_test_mappings"
kind = "donna.lib.run_script"
save_stdout_to = "mappings_stdout"
save_stderr_to = "mappings_stderr"
goto_on_success = "verify_governance_mappings"
goto_on_failure = "repair_source_test_mappings"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

relation_artifacts() {
    local relation="$1"
    local artifact="$2"
    depmesh -p llm dependencies --relation "$relation" "$artifact" |
        sed -n 's/^- \(@\/.*\)$/\1/p' |
        LC_ALL=C sort -u
}

assert_exact() {
    local relation="$1"
    local artifact="$2"
    shift 2
    local actual expected
    actual="$(relation_artifacts "$relation" "$artifact")"
    if (( $# )); then
        expected="$(printf '%s\n' "$@" | LC_ALL=C sort -u)"
    else
        expected=""
    fi
    if [[ "$actual" != "$expected" ]]; then
        printf 'Unexpected %s result for %s.\nExpected:\n%s\nActual:\n%s\n' \
            "$relation" "$artifact" "$expected" "$actual" >&2
        return 1
    fi
}

assert_contains() {
    local relation="$1"
    local artifact="$2"
    local expected="$3"
    if ! relation_artifacts "$relation" "$artifact" | rg -Fqx -- "$expected"; then
        printf '%s result for %s does not contain %s.\n' \
            "$relation" "$artifact" "$expected" >&2
        return 1
    fi
}

for source_path in R/*.R; do
    module="${source_path#R/}"
    test_path="tests/testthat/test-${module}"
    if [[ -f "$test_path" ]]; then
        assert_contains tested_by "@/${source_path}" "@/${test_path}"
        assert_contains tests "@/${test_path}" "@/${source_path}"
    fi
done

assert_exact tested_by @/R/normalize.R \
    @/tests/testthat/test-normalize.R
assert_exact tests @/tests/testthat/test-normalize.R \
    @/R/normalize.R

assert_exact tested_by @/R/correct_batch_effects.R \
    @/tests/testthat/test-batch_effect_steps.R \
    @/tests/testthat/test-correct_batch_effect.R
assert_exact tests @/tests/testthat/test-correct_batch_effect.R \
    @/R/batch_correction_helpers.R \
    @/R/correct_batch_effects.R

assert_exact tested_by @/R/correlation-based_diagnostics.R \
    @/tests/testthat/test-correlation_based_diagnostics.R
assert_exact tests @/tests/testthat/test-correlation_based_diagnostics.R \
    @/R/correlation-based_diagnostics.R

assert_exact tested_by @/R/pb_missing_filters.R \
    @/tests/testthat/test-pb_missing_helpers.R
assert_exact tests @/tests/testthat/test-pb_missing_helpers.R \
    @/R/ProBatchFeatures.R \
    @/R/pb_missing_filters.R

assert_exact tested_by @/R/ProBatchFeatures.R \
    @/tests/testthat/test-ProBatchFeatures.R \
    @/tests/testthat/test-batch_effect_steps.R \
    @/tests/testthat/test-pb_missing_helpers.R

assert_exact tested_by @/R/zzz_helpers.R \
    @/tests/testthat/test-batch_effect_steps.R
assert_exact tests @/tests/testthat/test-batch_effect_steps.R \
    @/R/ProBatchFeatures.R \
    @/R/correct_batch_effects.R \
    @/R/zzz_helpers.R

assert_exact tested_by @/R/proteome_wide_diagnostics.R \
    @/tests/testthat/test-explained_variance_plots.R \
    @/tests/testthat/test-proteome_wide_diagnostics.R
assert_exact tests @/tests/testthat/test-explained_variance_plots.R \
    @/R/proteome_wide_diagnostics.R

assert_exact tested_by @/R/pbf_input_helpers.R \
    @/tests/testthat/test-plot_helpers.R
assert_exact tests @/tests/testthat/test-plot_helpers.R \
    @/R/pbf_input_helpers.R \
    @/R/plot_helpers.R

assert_exact tested_by @/R/batch_correction_helpers.R \
    @/tests/testthat/test-correct_batch_effect.R

for source_path in R/*.R; do
    source_id="@/${source_path}"
    mapfile -t mapped_tests < <(relation_artifacts tested_by "$source_id")
    for test_id in "${mapped_tests[@]}"; do
        test_path="${test_id#@/}"
        [[ -f "$test_path" ]] || {
            printf 'Missing mapped test: %s\n' "$test_id" >&2
            exit 1
        }
        assert_contains tests "$test_id" "$source_id"
    done
done

for test_path in tests/testthat/test-*.R; do
    test_id="@/${test_path}"
    mapfile -t mapped_sources < <(relation_artifacts tests "$test_id")
    for source_id in "${mapped_sources[@]}"; do
        source_path="${source_id#@/}"
        [[ -f "$source_path" ]] || {
            printf 'Missing mapped source: %s\n' "$source_id" >&2
            exit 1
        }
        assert_contains tested_by "$source_id" "$test_id"
    done
done

printf 'Validated conventional mappings, explicit filename exceptions, integration mappings, and reciprocal existing paths.\n'
```

## Repair Source and Test Mappings

```toml donna
id = "repair_source_test_mappings"
kind = "donna.lib.request_action"
```

The source/test relation gate failed.

Standard output:

```text
{{ donna.lib.task_variable("mappings_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("mappings_stderr") }}
```

1. Inspect the named R sources and tests to confirm direct or intentional integration ownership.
2. Query all relations for each artifact before editing `{{ donna.lib.path("@/depmesh.toml") }}` or a helper.
3. Keep conventional rules scalable and declare every exception in both directions.
4. Repair only the relation model or an incorrect verifier expectation. Do not change package behavior or inspect `man/`.
5. After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`.
6. If ownership is ambiguous or a specification conflict remains, `{{ donna.lib.goto("blocked") }}`.

## Verify Governance Mappings

```toml donna
id = "verify_governance_mappings"
kind = "donna.lib.run_script"
save_stdout_to = "governance_stdout"
save_stderr_to = "governance_stderr"
goto_on_success = "verify_path_safety"
goto_on_failure = "repair_governance_mappings"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
shopt -s globstar nullglob

relation_artifacts() {
    local relation="$1"
    local artifact="$2"
    depmesh -p llm dependencies --relation "$relation" "$artifact" |
        sed -n 's/^- \(@\/.*\)$/\1/p' |
        LC_ALL=C sort -u
}

assert_exact() {
    local relation="$1"
    local artifact="$2"
    shift 2
    local actual expected
    actual="$(relation_artifacts "$relation" "$artifact")"
    if (( $# )); then
        expected="$(printf '%s\n' "$@" | LC_ALL=C sort -u)"
    else
        expected=""
    fi
    if [[ "$actual" != "$expected" ]]; then
        printf 'Unexpected %s result for %s.\nExpected:\n%s\nActual:\n%s\n' \
            "$relation" "$artifact" "$expected" "$actual" >&2
        return 1
    fi
}

assert_contains() {
    local relation="$1"
    local artifact="$2"
    local expected="$3"
    if ! relation_artifacts "$relation" "$artifact" | rg -Fqx -- "$expected"; then
        printf '%s result for %s does not contain %s.\n' \
            "$relation" "$artifact" "$expected" >&2
        return 1
    fi
}

spec_paths=(specs/**/*.md)
spec_ids=()
meta_governed_ids=("@/AGENTS.md")
for path in "${spec_paths[@]}"; do
    id="@/${path}"
    spec_ids+=("$id")
    if [[ "$path" != "specs/meta/general.md" ]]; then
        meta_governed_ids+=("$id")
        assert_contains governed_by "$id" @/specs/meta/general.md
    else
        assert_exact governed_by "$id"
    fi
    assert_contains governs "$id" @/AGENTS.md
done

assert_exact governed_by @/AGENTS.md "${spec_ids[@]}"
assert_exact governs @/specs/meta/general.md "${meta_governed_ids[@]}"

helper_paths=()
for path in bin/depemesh/**/*; do
    [[ -f "$path" ]] && helper_paths+=("$path")
done

behavior_paths=(
    .Rbuildignore
    AGENTS.md
    depmesh.toml
    donna.toml
    R/*.R
    tests/**/*.R
    vignettes/*.Rmd
    workflows/*.donna.md
    workflows/reports/*.md
    "${helper_paths[@]}"
)
behavior_ids=()
for path in "${behavior_paths[@]}"; do
    behavior_ids+=("@/${path}")
    assert_contains governed_by "@/${path}" @/specs/behavior/files_relations.md
done
assert_exact governs @/specs/behavior/files_relations.md "${behavior_ids[@]}"

workflow_paths=(
    AGENTS.md
    donna.toml
    workflows/*.donna.md
    workflows/reports/*.md
)
workflow_ids=()
for path in "${workflow_paths[@]}"; do
    workflow_ids+=("@/${path}")
    assert_contains governed_by "@/${path}" @/specs/general/workflows.md
done
assert_exact governs @/specs/general/workflows.md "${workflow_ids[@]}"

assert_exact governed_by @/.session/donna/state.json

printf 'Validated meta, relation-contract, workflow, helper, and session-state governance in both directions.\n'
```

## Repair Governance Mappings

```toml donna
id = "repair_governance_mappings"
kind = "donna.lib.request_action"
```

The governance relation gate failed.

Standard output:

```text
{{ donna.lib.task_variable("governance_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("governance_stderr") }}
```

Query the affected artifacts in both directions, then repair only the scalable governance rules, indexed specification paths, workflow/catalog synchronization, or verifier expectation. Keep session artifacts absent from permanent dependencies and do not inspect `man/`.

After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`. If the required ownership is unclear or conflicts with a specification, `{{ donna.lib.goto("blocked") }}`.

## Verify Path Safety

```toml donna
id = "verify_path_safety"
kind = "donna.lib.run_script"
save_stdout_to = "safety_stdout"
save_stderr_to = "safety_stderr"
goto_on_success = "verify_workflow_artifact"
goto_on_failure = "repair_path_safety"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
shopt -s globstar nullglob

project_root="$(pwd -P)"
helper_paths=()
for path in bin/depemesh/**/*; do
    [[ -f "$path" ]] && helper_paths+=("$path")
done

input_paths=(
    .Rbuildignore
    AGENTS.md
    depmesh.toml
    donna.toml
    specs/**/*.md
    R/*.R
    tests/**/*.R
    vignettes/*.Rmd
    workflows/*.donna.md
    workflows/reports/*.md
    "${helper_paths[@]}"
    .session/donna/state.json
)

input_ids=()
for path in "${input_paths[@]}"; do
    input_ids+=("@/${path}")
done
mapfile -t input_ids < <(printf '%s\n' "${input_ids[@]}" | LC_ALL=C sort -u)

edge_count=0
for input_id in "${input_ids[@]}"; do
    query_output="$(depmesh -p llm dependencies "$input_id")"
    while IFS= read -r dependency; do
        [[ -n "$dependency" ]] || continue
        [[ "$dependency" == @/* ]] || {
            printf 'Non-root-anchored dependency returned for %s: %s\n' \
                "$input_id" "$dependency" >&2
            exit 1
        }
        [[ "$dependency" != "@/man" && "$dependency" != "@/man/"* ]] || {
            printf 'Generated manual dependency returned for %s: %s\n' \
                "$input_id" "$dependency" >&2
            exit 1
        }
        [[ "$dependency" != *"/../"* && "$dependency" != *"/.." ]] || {
            printf 'Parent traversal returned for %s: %s\n' \
                "$input_id" "$dependency" >&2
            exit 1
        }
        relative="${dependency#@/}"
        [[ -e "$relative" ]] || {
            printf 'Missing dependency returned for %s: %s\n' \
                "$input_id" "$dependency" >&2
            exit 1
        }
        resolved="$(realpath -e -- "$relative")"
        case "$resolved" in
            "$project_root"/*) ;;
            *)
                printf 'Dependency resolves outside the project for %s: %s\n' \
                    "$input_id" "$dependency" >&2
                exit 1
                ;;
        esac
        edge_count=$((edge_count + 1))
    done < <(printf '%s\n' "$query_output" | sed -n 's/^- //p')
done

relation_sources=(depmesh.toml "${helper_paths[@]}")
if rg -n '(^|[^[:alnum:]_])man(/|[^[:alnum:]_]|$)' "${relation_sources[@]}"; then
    printf 'A relation rule or helper contains a generated-manual path token.\n' >&2
    exit 1
fi

required_build_exclusions=(
    '^AGENTS\.md$'
    '^bin$'
    '^depmesh\.toml$'
    '^donna\.toml$'
    '^specs$'
    '^workflows$'
    '^\.session$'
)
for pattern in "${required_build_exclusions[@]}"; do
    if ! rg -Fqx -- "$pattern" .Rbuildignore; then
        printf 'Missing .Rbuildignore entry: %s\n' "$pattern" >&2
        exit 1
    fi
done

input_count="$(printf '%s\n' "${input_ids[@]}" | wc -l | tr -d ' ')"
printf 'Validated %d allowed inputs and %d project-local dependency edges without accessing generated manuals.\n' \
    "$input_count" "$edge_count"
```

## Repair Path Safety

```toml donna
id = "repair_path_safety"
kind = "donna.lib.request_action"
```

The project-local path and generated-manual exclusion gate failed.

Standard output:

```text
{{ donna.lib.task_variable("safety_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("safety_stderr") }}
```

Inspect only `{{ donna.lib.path("@/depmesh.toml") }}`, helper sources under `{{ donna.lib.path("@/bin/depemesh") }}`, and the named project artifact. Do not inspect the forbidden destination. Repair broad patterns, helper traversal, non-project outputs, missing outputs, or build exclusions without weakening the relation contract.

After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`. If safe repair requires a new policy decision, `{{ donna.lib.goto("blocked") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "review_relation_sources"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/verify-file-relations.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"
rg -Fq -- '- Artifact: `./workflows/verify-file-relations.donna.md`' \
    specs/general/workflows.md

catalog_state="$(
    awk '
        /^### Verify file relations$/ { in_section = 1; next }
        in_section && /^### / { exit }
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
        printf 'Unexpected workflow catalog state: %s\n' "$catalog_state" >&2
        exit 1
        ;;
esac

depmesh -p llm dependencies "$workflow_id"

git diff --check -- \
    AGENTS.md \
    depmesh.toml \
    specs/general/workflows.md \
    "$workflow_path"

if rg -n '[[:blank:]]+$' "$workflow_path"; then
    printf 'Workflow contains trailing whitespace.\n' >&2
    exit 1
fi
if [[ -n "$(tail -c 1 "$workflow_path")" ]]; then
    printf 'Workflow must end with a newline.\n' >&2
    exit 1
fi

printf 'Workflow source is governed and synchronized with its catalog entry; validation and rendering remain agent-side checks outside an executing Donna process.\n'
```

## Repair Workflow Artifact

```toml donna
id = "repair_workflow_artifact"
kind = "donna.lib.request_action"
```

The workflow-artifact gate failed.

Standard output:

```text
{{ donna.lib.task_variable("workflow_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("workflow_stderr") }}
```

Repair this workflow, its catalog entry, affected agent guidance, or its scalable governance relation. Keep the catalog state at `implementation in progress` until the first representative execution finishes. Validate any control-flow change before continuing.

After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`. If repair requires developer direction, `{{ donna.lib.goto("blocked") }}`.

## Review Relation Sources

```toml donna
id = "review_relation_sources"
kind = "donna.lib.request_action"
```

Deterministic relation gates passed.

1. Inspect the exact in-scope diff and untracked workflow without modifying the Git index.
2. Read `{{ donna.lib.path("@/depmesh.toml") }}` and every helper under `{{ donna.lib.path("@/bin/depemesh") }}` completely. Do not inspect `man/`.
3. Confirm relation definitions remain directional and return-oriented; conventional rules precede explicit exceptions; every intentional exception is bidirectional; and reverse results match forward ownership.
4. Confirm every input and output pattern is project-local. Ensure no broad pattern or helper can read, enumerate, return, or check generated manuals, even when no literal forbidden path appears.
5. Confirm helpers, if any, are deterministic, read-only, non-networked, independent of user state, executable from the project root, and emit only the specified protocol.
6. While Donna is awaiting this action and no workflow operation is executing, run `donna -p llm validate @/workflows/verify-file-relations.donna.md` and `donna -p llm render @/workflows/verify-file-relations.donna.md --mode view`; review the rendered instructions.
7. Confirm the catalog, agent guidance, Depmesh rules, and this workflow form one coherent implementation change.
8. If the complete relation contract is satisfied, `{{ donna.lib.goto("finish") }}`.
9. If a repairable discrepancy remains, `{{ donna.lib.goto("repair_review_findings") }}`.
10. If the specifications conflict or ownership needs a developer decision, `{{ donna.lib.goto("blocked") }}`.

## Repair Review Findings

```toml donna
id = "repair_review_findings"
kind = "donna.lib.request_action"
```

Query all relations for every affected artifact, then repair the reviewed rule, helper, workflow, catalog, or guidance discrepancy. Preserve unrelated changes, keep the Git index untouched, and do not inspect generated manuals.

After repair, `{{ donna.lib.goto("verify_relation_definitions") }}`. If the finding cannot be resolved within the specified relation contract, `{{ donna.lib.goto("blocked") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The complete file-relation contract passed deterministic checks and source-level review. Report relation definitions, bidirectional source/test and governance coverage, path-safety evidence, files repaired, and any intentional empty mappings.

If this was the first representative execution and the catalog still says `implementation in progress`, mark it `implemented` only after this finish, then rerun workflow validation, view rendering, governance queries, and scoped diff checks before reporting completion.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The file-relation workflow cannot continue without developer input. Report completed gates, the exact relation or specification conflict, affected artifacts, and the decision required. Do not inspect generated manuals, stage files, create commits, or reset the Donna session.
