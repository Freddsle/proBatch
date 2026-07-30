# Review Documentation Sources

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Review public API documentation at its editable Roxygen2 source under `R/`, verify it against the affected interfaces, and hand generated documentation and namespace regeneration to the maintainer.

This workflow MUST NOT access `man/` or `NAMESPACE`, generate documentation, run package-wide checks, use the network, change executable package behavior, stage files, create commits, or reset the Donna session.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/meta/general.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}` completely.
2. Read the built-in Donna workflow and Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations separately for this workflow, its catalog, `{{ donna.lib.path("@/AGENTS.md") }}`, and every R source already known to be in the initiating scope.
3. Inspect the non-generated worktree without modifying the Git index. Preserve unrelated user changes and distinguish initiating changes from pre-existing changes.
4. Confirm this action belongs to the current Donna task and no unrelated action request is being displaced. Do not reset the session.
5. Confirm this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna operation was executing.
6. Confirm the authorized work is limited to source-level documentation review and Roxygen2 comment repair. Do not modify executable R code, ordinary comments, generated output, package behavior, or files outside the reviewed scope.
7. Do not read, enumerate, compare, hash, validate, or otherwise inspect `man/` or `NAMESPACE`. Do not run documentation generation, spelling checks that inspect package manuals, package installation, package building, or package-wide checks.
8. If the governed scope and source-only boundary are understood, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
9. If specifications conflict, the initiating scope is unavailable, or the requested work requires generated output, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "discover_documentation_scope"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/review-documentation-sources.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"

artifact_entry_count="$(
    rg -Fxc -- '- Artifact: `./workflows/review-documentation-sources.donna.md`' \
        specs/general/workflows.md
)"
if (( artifact_entry_count != 1 )); then
    printf 'Workflow catalog must contain exactly one documentation-source artifact entry; found %d.\n' \
        "$artifact_entry_count" >&2
    exit 1
fi

catalog_state="$(
    awk '
        /^### Review documentation sources$/ {
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
        printf 'Unexpected documentation-source workflow catalog state: %s\n' \
            "$catalog_state" >&2
        exit 1
        ;;
esac

if ! rg -Fq -- '@/workflows/review-documentation-sources.donna.md' AGENTS.md; then
    printf 'AGENTS.md must route documentation-source triggers to this workflow.\n' >&2
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

printf 'Workflow source, lifecycle state, agent routing, package exclusion, governance, and non-generated diff hygiene are synchronized.\n'
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

## Discover Documentation Scope

```toml donna
id = "discover_documentation_scope"
kind = "donna.lib.run_script"
save_stdout_to = "discovery_stdout"
save_stderr_to = "discovery_stderr"
goto_on_success = "review_documentation_scope"
goto_on_failure = "repair_scope_discovery"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/review-documentation-sources"
scope_file="$plan_dir/scope.txt"
deleted_file="$plan_dir/deleted-sources.txt"
plan_file="$plan_dir/plan.md"
generation_pointer="$plan_dir/generation-current"

if [[ -L "$plan_dir" ]]; then
    printf 'Session planning directory must not be a symlink.\n' >&2
    exit 1
fi
install -d -m 700 -- "$plan_dir"

for output_path in \
    "$scope_file" \
    "$deleted_file" \
    "$plan_file" \
    "$generation_pointer"
do
    if [[ -L "$output_path" ]]; then
        printf 'Session planning output must not be a symlink: %s\n' \
            "$output_path" >&2
        exit 1
    fi
done

scope_temp="$(mktemp "$plan_dir/scope.XXXXXXXX")"
deleted_temp="$(mktemp "$plan_dir/deleted-sources.XXXXXXXX")"
plan_temp="$(mktemp "$plan_dir/plan.XXXXXXXX")"
generation_path="$(mktemp "$plan_dir/generation.XXXXXXXX")"
generation_temp="$(mktemp "$plan_dir/generation-current.XXXXXXXX")"
trap 'rm -f -- "$scope_temp" "$deleted_temp" "$plan_temp" "$generation_temp"' EXIT

{
    git --no-optional-locks diff --name-only --diff-filter=ACMRTUXB -- R
    git --no-optional-locks diff --cached --name-only \
        --diff-filter=ACMRTUXB -- R
    git --no-optional-locks ls-files --others --exclude-standard -- R
} |
    awk '/^R\/[^/]+[.]R$/ { print "@/" $0 }' |
    LC_ALL=C sort -u > "$scope_temp"

{
    git --no-optional-locks diff --name-only --diff-filter=D -- R
    git --no-optional-locks diff --cached --name-only --diff-filter=D -- R
} |
    awk '/^R\/[^/]+[.]R$/ { print "@/" $0 }' |
    LC_ALL=C sort -u > "$deleted_temp"

candidate_count="$(wc -l < "$scope_temp" | tr -d ' ')"
deleted_count="$(wc -l < "$deleted_temp" | tr -d ' ')"

{
    printf '# Documentation-source review plan\n\n'
    printf '## Scope rationale\n\n'
    if (( candidate_count )); then
        printf 'Candidate R sources were derived from current staged, unstaged, and non-ignored untracked worktree paths. They are candidates only and must be pruned to the initiating documentation-review scope.\n\n'
    else
        printf 'No changed R source was discoverable. Add only an explicitly requested source or a deliberately selected representative source before validation.\n\n'
    fi
    printf '## Review notes\n\n'
    printf -- '- Confirm every selected source owns editable Roxygen2 documentation for an affected public interface or documented dataset.\n'
    printf -- '- Record why every candidate is included and why any worktree candidate is removed.\n'
    if (( deleted_count )); then
        printf -- '- Deleted R sources were detected and are recorded separately. A deleted source cannot enter the editable-source scope; determine whether the deletion requires developer direction.\n'
    fi
} > "$plan_temp"

printf '%s\n' "${generation_path##*/}" > "$generation_temp"
mv -f -- "$scope_temp" "$scope_file"
mv -f -- "$deleted_temp" "$deleted_file"
mv -f -- "$plan_temp" "$plan_file"
mv -f -- "$generation_temp" "$generation_pointer"

printf 'CANDIDATE_SCOPE_COUNT=%s\nDELETED_SOURCE_COUNT=%s\nGENERATION=%s\n' \
    "$candidate_count" "$deleted_count" "${generation_path##*/}"
```

## Repair Scope Discovery

```toml donna
id = "repair_scope_discovery"
kind = "donna.lib.request_action"
```

Documentation-scope discovery failed.

Standard output:

```text
{{ donna.lib.task_variable("discovery_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("discovery_stderr") }}
```

1. Diagnose only the discovery command, repository state, session directory, or an incorrect path expectation.
2. Preserve project artifacts, unrelated changes, and the Git index.
3. If the discovery implementation can be repaired in this workflow, query its relations, repair it, validate the control flow, and `{{ donna.lib.goto("verify_workflow_artifact") }}`.
4. If safe discovery requires generated output, a repository-state decision, or developer input, `{{ donna.lib.goto("record_blocker") }}`.

## Review Documentation Scope

```toml donna
id = "review_documentation_scope"
kind = "donna.lib.request_action"
```

Candidate discovery completed.

```text
{{ donna.lib.task_variable("discovery_stdout") }}
```

1. Inspect `{{ donna.lib.path("@/.session/donna/review-documentation-sources/scope.txt") }}`, `{{ donna.lib.path("@/.session/donna/review-documentation-sources/deleted-sources.txt") }}`, and `{{ donna.lib.path("@/.session/donna/review-documentation-sources/plan.md") }}`.
2. Treat worktree-derived entries as candidates, not authorization. Remove unrelated or pre-existing sources and add initiating sources known from the request but absent from discovery.
3. For a first representative execution with no initiating R change, select one small existing public-interface source and state why it exercises the workflow meaningfully.
4. Before adding a source, run `depmesh -p llm dependencies` for its root-anchored artifact identifier and inspect every returned relation. Use mapped tests only as read-only behavioral evidence.
5. Keep one sorted, unique `@/R/<file>.R` identifier per line. Select existing regular files only; never add `NAMESPACE`, generated output, directories, deleted sources, tests, vignettes, or files outside `R/`.
6. Update the plan with the exact scope rationale, interface or dataset anchors to review, relation evidence, and known exclusions. At least one source is required.
7. Edit only these session-local planning files during this action.
8. If the reviewed scope is authorized and complete, `{{ donna.lib.goto("validate_documentation_scope") }}`.
9. If the planning files need another review, `{{ donna.lib.goto("review_documentation_scope") }}`.
10. If a deleted interface, conflicting requirement, or missing authority prevents a safe source review, `{{ donna.lib.goto("record_blocker") }}`.

## Validate Documentation Scope

```toml donna
id = "validate_documentation_scope"
kind = "donna.lib.run_script"
save_stdout_to = "scope_stdout"
save_stderr_to = "scope_stderr"
goto_on_success = "review_documentation_sources"
goto_on_failure = "repair_documentation_scope"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/review-documentation-sources"
scope_file="$plan_dir/scope.txt"
plan_file="$plan_dir/plan.md"
generation_pointer="$plan_dir/generation-current"
baseline_pointer="$plan_dir/baseline-current"

test -s "$scope_file"
test -s "$plan_file"
test -s "$generation_pointer"
for input_path in "$plan_dir" "$scope_file" "$plan_file" "$generation_pointer"; do
    if [[ -L "$input_path" ]]; then
        printf 'Session planning input must not be a symlink: %s\n' \
            "$input_path" >&2
        exit 1
    fi
done
rg -Fq -- '## Scope rationale' "$plan_file"

generation_name="$(sed -n '1p' "$generation_pointer")"
generation_path="$plan_dir/$generation_name"
if [[ ! "$generation_name" =~ ^generation[.][A-Za-z0-9]+$ ]] ||
   [[ ! -f "$generation_path" ]] || [[ -L "$generation_path" ]]; then
    printf 'Current discovery generation is missing or malformed.\n' >&2
    exit 1
fi

sorted_scope="$(LC_ALL=C sort -u -- "$scope_file")"
if [[ "$sorted_scope" != "$(sed -n '/./p' "$scope_file")" ]]; then
    printf 'Scope identifiers must be nonblank, sorted, and unique.\n' >&2
    exit 1
fi

relation_output_file="$(mktemp "$plan_dir/relations.XXXXXXXX")"
trap 'rm -f -- "$relation_output_file"' EXIT

scope_paths=()
while IFS= read -r artifact_id; do
    if [[ ! "$artifact_id" =~ ^@/R/[^/[:space:]]+[.]R$ ]]; then
        printf 'Invalid documentation-source scope identifier: %s\n' \
            "$artifact_id" >&2
        exit 1
    fi
    source_path="${artifact_id#@/}"
    if [[ ! -f "$source_path" || -L "$source_path" ]]; then
        printf 'Scoped source must be an existing regular non-symlink file: %s\n' \
            "$source_path" >&2
        exit 1
    fi
    scope_paths+=("$source_path")

    depmesh -p llm dependencies "$artifact_id" > "$relation_output_file"
    if rg -Fqx -- '## warnings' "$relation_output_file"; then
        printf 'Depmesh returned warnings for %s:\n' "$artifact_id" >&2
        sed -n '1,240p' "$relation_output_file" >&2
        exit 1
    fi
    if sed -n 's/^- \(@\/.*\)$/\1/p' "$relation_output_file" |
        rg -q '^@/(man(?:/|$)|NAMESPACE$)'; then
        printf 'Dependency discovery returned forbidden generated output for %s.\n' \
            "$artifact_id" >&2
        exit 1
    fi
done < "$scope_file"

scope_count="$(wc -l < "$scope_file" | tr -d ' ')"
if (( scope_count == 0 )); then
    printf 'At least one existing R source is required.\n' >&2
    exit 1
fi

emit_worktree_manifest() {
    git --no-optional-locks ls-files -z --cached --others --exclude-standard -- "$@" |
        LC_ALL=C sort -z |
        while IFS= read -r -d '' path; do
            printf '%s\0' "$path" || exit $?
            if [[ -L "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'symlink\0%s\0' "$mode" || exit $?
                readlink -z -- "$path" |
                    sha256sum |
                    awk '{ printf "%s%c", $1, 0 }'
            elif [[ -f "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'file\0%s\0' "$mode" || exit $?
                sha256sum < "$path" |
                    awk '{ printf "%s%c", $1, 0 }'
            elif [[ -e "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'other\0%s\0-\0' "$mode" || exit $?
            else
                printf 'missing\0-\0-\0' || exit $?
            fi
        done
}

emit_index_manifest() {
    printf 'index-stage-v1\0'
    git --no-optional-locks ls-files --stage -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-assume-skip-v1\0'
    git --no-optional-locks ls-files --stage -v -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-fsmonitor-v1\0'
    git --no-optional-locks ls-files --stage -f -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-entry-flags-v1\0'
    LC_ALL=C git --no-optional-locks ls-files --stage --debug -z -- "$@"
    printf '\0index-resolve-undo-v1\0'
    git --no-optional-locks ls-files --resolve-undo -z -- "$@" |
        LC_ALL=C sort -z
}

if [[ -e "$baseline_pointer" ]]; then
    if [[ -L "$baseline_pointer" || ! -f "$baseline_pointer" ]]; then
        printf 'Current baseline pointer must be a regular non-symlink file.\n' >&2
        exit 1
    fi
    previous_attempt_name="$(sed -n '1p' "$baseline_pointer")"
    previous_attempt_dir="$plan_dir/$previous_attempt_name"
    if [[ ! "$previous_attempt_name" =~ ^baseline[.][A-Za-z0-9]+$ ]] ||
       [[ ! -d "$previous_attempt_dir" ]] ||
       [[ -L "$previous_attempt_dir" ]] ||
       [[ ! -f "$previous_attempt_dir/state" ]] ||
       [[ -L "$previous_attempt_dir/state" ]] ||
       [[ ! -f "$previous_attempt_dir/source-baseline.tsv" ]] ||
       [[ -L "$previous_attempt_dir/source-baseline.tsv" ]]; then
        printf 'Existing baseline attempt is missing or malformed.\n' >&2
        exit 1
    fi
    previous_generation="$(
        sed -n 's/^generation=//p' "$previous_attempt_dir/state"
    )"
    if [[ "$previous_generation" == "$generation_name" ]]; then
        while IFS=$'\t' read -r previous_source _; do
            previous_copy="$previous_attempt_dir/$previous_source"
            if [[ ! -f "$previous_copy" ]] || [[ -L "$previous_copy" ]] ||
               [[ ! -f "$previous_source" ]] || [[ -L "$previous_source" ]] ||
               ! cmp -s -- "$previous_copy" "$previous_source"; then
                printf 'Scope cannot be re-baselined after a scoped source changed: %s\n' \
                    "$previous_source" >&2
                exit 1
            fi
        done < "$previous_attempt_dir/source-baseline.tsv"
    fi
fi

attempt_dir="$(mktemp -d "$plan_dir/baseline.XXXXXXXX")"
baseline_manifest="$attempt_dir/source-baseline.tsv"

for source_path in "${scope_paths[@]}"; do
    target_dir="$attempt_dir/${source_path%/*}"
    install -d -m 700 -- "$target_dir"
    cp -p -- "$source_path" "$attempt_dir/$source_path"
    mode="$(stat -c '%f' -- "$source_path")"
    projection_sha="$(
        awk '!/^#'\''/' "$source_path" |
            sha256sum |
            awk '{ print $1 }'
    )"
    printf '%s\t%s\t%s\n' "$source_path" "$mode" "$projection_sha" \
        >> "$baseline_manifest"
done

protected_pathspecs=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
for source_path in "${scope_paths[@]}"; do
    protected_pathspecs+=(":(top,exclude)${source_path}")
done

index_fingerprint="$(
    emit_index_manifest \
        . \
        ':(top,exclude)man/**' \
        ':(top,exclude)NAMESPACE' |
        sha256sum |
        awk '{ print $1 }'
)"
protected_fingerprint="$(
    {
        git --no-optional-locks status --porcelain=v2 -z \
            --untracked-files=all -- "${protected_pathspecs[@]}" &&
            printf '\0protected-worktree-v1\0' &&
            emit_worktree_manifest "${protected_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"
head_revision="$(git --no-optional-locks rev-parse --verify 'HEAD^{commit}')"
scope_sha="$(sha256sum < "$scope_file" | awk '{ print $1 }')"

if [[ ! "$head_revision" =~ ^[0-9a-f]{40,64}$ ]] ||
   [[ ! "$index_fingerprint" =~ ^[0-9a-f]{64}$ ]] ||
   [[ ! "$protected_fingerprint" =~ ^[0-9a-f]{64}$ ]] ||
   [[ ! "$scope_sha" =~ ^[0-9a-f]{64}$ ]]; then
    printf 'Scope validation produced malformed repository fingerprints.\n' >&2
    exit 1
fi

{
    printf 'version=1\n'
    printf 'head=%s\n' "$head_revision"
    printf 'index=%s\n' "$index_fingerprint"
    printf 'protected=%s\n' "$protected_fingerprint"
    printf 'scope=%s\n' "$scope_sha"
    printf 'generation=%s\n' "$generation_name"
} > "$attempt_dir/state"
baseline_pointer_temp="$(mktemp "$plan_dir/baseline-current.XXXXXXXX")"
printf '%s\n' "${attempt_dir##*/}" > "$baseline_pointer_temp"
mv -f -- "$baseline_pointer_temp" "$baseline_pointer"

printf 'SCOPE_COUNT=%d\nBASELINE_ATTEMPT=%s\nHEAD=%s\nINDEX_SHA256=%s\nPROTECTED_SHA256=%s\nSCOPE_SHA256=%s\n' \
    "$scope_count" "${attempt_dir##*/}" "$head_revision" \
    "$index_fingerprint" "$protected_fingerprint" "$scope_sha"
```

## Repair Documentation Scope

```toml donna
id = "repair_documentation_scope"
kind = "donna.lib.request_action"
```

Documentation-scope validation failed.

Standard output:

```text
{{ donna.lib.task_variable("scope_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("scope_stderr") }}
```

1. If a session-local identifier, rationale, ordering, or stale candidate is wrong and every source in the active baseline remains byte-identical to its saved copy, edit only the planning files and `{{ donna.lib.goto("validate_documentation_scope") }}`.
2. Query all relations before adding or replacing a source.
3. If the workflow verifier or scalable governance mapping is wrong, repair it, validate the control flow, and `{{ donna.lib.goto("verify_workflow_artifact") }}`.
4. Preserve project artifacts, unrelated changes, generated output, and the Git index.
5. If re-baselining would conceal a scoped source edit, or a correct nonempty scope cannot be established safely, `{{ donna.lib.goto("record_blocker") }}`.

## Review Documentation Sources

```toml donna
id = "review_documentation_sources"
kind = "donna.lib.request_action"
```

Validated documentation scope and protected baseline:

```text
{{ donna.lib.task_variable("scope_stdout") }}
```

1. Read every scoped R source completely and inspect the contiguous Roxygen2 block for each affected exported function, class, S3 or S4 method, generic, or documented dataset. Do not inspect generated output.
2. Query all configured Depmesh relations separately for every scoped source. Read mapped tests only as behavioral evidence; do not edit or execute them in this workflow.
3. Compare documentation with current formals and defaults, accepted types and shapes, return value and side effects, errors and deprecations, examples, shared-topic links, inherited parameters, export directives, and source-level S3 or S4 registration as applicable.
4. For datasets, review the editable title, description, usage, format, source, and field descriptions against the packaged data contract available outside generated output.
5. Improve only `#'` Roxygen2 comment lines in the authorized sources. Do not change executable R code, ordinary comments, whitespace outside Roxygen blocks, tests, vignettes, package metadata, `NAMESPACE`, or any generated file.
6. Do not run Roxygen2, devtools documentation commands, spelling checks over the package, examples, package installation, package building, or package-wide checks.
7. If the source-level documentation is ready for deterministic verification, `{{ donna.lib.goto("verify_documentation_sources") }}`.
8. If the scope must change before any scoped source differs from its saved baseline, update only the session plan and `{{ donna.lib.goto("validate_documentation_scope") }}` to capture a new baseline before further edits.
9. If the scope must change after a scoped source differs from its baseline, preserve the edit and `{{ donna.lib.goto("record_blocker") }}` rather than replacing its evidence boundary.
10. If correct documentation requires a behavior change, generated-output inspection, or a product decision, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Documentation Sources

```toml donna
id = "verify_documentation_sources"
kind = "donna.lib.run_script"
save_stdout_to = "source_stdout"
save_stderr_to = "source_stderr"
goto_on_success = "review_documentation_semantics"
goto_on_failure = "repair_documentation_sources"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/review-documentation-sources"
scope_file="$plan_dir/scope.txt"
baseline_pointer="$plan_dir/baseline-current"

if [[ ! -f "$scope_file" ]] || [[ -L "$scope_file" ]] ||
   [[ ! -f "$baseline_pointer" ]] || [[ -L "$baseline_pointer" ]]; then
    printf 'Scope or baseline pointer is missing or unsafe.\n' >&2
    exit 1
fi

attempt_name="$(sed -n '1p' "$baseline_pointer")"
attempt_dir="$plan_dir/$attempt_name"

if [[ ! "$attempt_name" =~ ^baseline[.][A-Za-z0-9]+$ ]] ||
   [[ ! -d "$attempt_dir" ]] || [[ -L "$attempt_dir" ]] ||
   [[ ! -f "$attempt_dir/state" ]] || [[ -L "$attempt_dir/state" ]] ||
   [[ ! -f "$attempt_dir/source-baseline.tsv" ]] ||
   [[ -L "$attempt_dir/source-baseline.tsv" ]]; then
    printf 'Current baseline attempt is missing or malformed.\n' >&2
    exit 1
fi

state_value() {
    local key="$1"
    sed -n "s/^${key}=//p" "$attempt_dir/state"
}

scope_sha="$(sha256sum < "$scope_file" | awk '{ print $1 }')"
if [[ "$scope_sha" != "$(state_value scope)" ]]; then
    printf 'Documentation scope changed after baseline capture.\n' >&2
    exit 1
fi

scope_paths=()
while IFS= read -r artifact_id; do
    source_path="${artifact_id#@/}"
    scope_paths+=("$source_path")

    baseline_path="$attempt_dir/$source_path"
    if [[ ! -f "$baseline_path" ]] || [[ -L "$baseline_path" ]] ||
       [[ ! -f "$source_path" ]] || [[ -L "$source_path" ]]; then
        printf 'Scoped source or saved copy is missing or unsafe: %s\n' \
            "$source_path" >&2
        exit 1
    fi

    expected_record="$(
        awk -F '\t' -v path="$source_path" '$1 == path { print; exit }' \
            "$attempt_dir/source-baseline.tsv"
    )"
    if [[ -z "$expected_record" ]]; then
        printf 'No source baseline record exists for %s.\n' "$source_path" >&2
        exit 1
    fi
    IFS=$'\t' read -r _ expected_mode expected_projection \
        <<< "$expected_record"
    actual_mode="$(stat -c '%f' -- "$source_path")"
    actual_projection="$(
        awk '!/^#'\''/' "$source_path" |
            sha256sum |
            awk '{ print $1 }'
    )"
    if [[ "$actual_mode" != "$expected_mode" ]] ||
       [[ "$actual_projection" != "$expected_projection" ]]; then
        printf 'Non-Roxygen source content or file mode changed for %s.\n' \
            "$source_path" >&2
        exit 1
    fi

    diff_output="$(mktemp "$plan_dir/source-diff.XXXXXXXX")"
    diff_status=0
    diff -u -- "$baseline_path" "$source_path" > "$diff_output" ||
        diff_status=$?
    if (( diff_status > 1 )); then
        printf 'Could not compare source baseline for %s.\n' \
            "$source_path" >&2
        rm -f -- "$diff_output"
        exit "$diff_status"
    fi
    while IFS= read -r diff_line; do
        case "$diff_line" in
            '--- '*|'+++ '*|'@@ '*) ;;
            +*|-*)
                changed_line="${diff_line:1}"
                if [[ "$changed_line" != "#'"* ]]; then
                    printf 'Documentation repair changed a non-Roxygen line in %s: %s\n' \
                        "$source_path" "$diff_line" >&2
                    rm -f -- "$diff_output"
                    exit 1
                fi
                ;;
        esac
    done < "$diff_output"
    rm -f -- "$diff_output"

    relation_output="$(depmesh -p llm dependencies "$artifact_id")"
    if printf '%s\n' "$relation_output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for %s:\n%s\n' \
            "$artifact_id" "$relation_output" >&2
        exit 1
    fi
    if printf '%s\n' "$relation_output" |
        sed -n 's/^- \(@\/.*\)$/\1/p' |
        rg -q '^@/(man(?:/|$)|NAMESPACE$)'; then
        printf 'Dependency discovery returned forbidden generated output for %s.\n' \
            "$artifact_id" >&2
        exit 1
    fi
done < "$scope_file"

DOC_REVIEW_SCOPE_FILE="$scope_file" \
    pixi run --as-is Rscript --vanilla - <<'RSCRIPT'
scope_file <- Sys.getenv("DOC_REVIEW_SCOPE_FILE")
scope_ids <- readLines(scope_file, warn = FALSE)
for (artifact_id in scope_ids) {
    source_path <- sub("^@/", "", artifact_id)
    invisible(parse(file = source_path, keep.source = TRUE))
    lines <- readLines(source_path, warn = FALSE)
    roxygen_lines <- grep("^#'", lines, value = TRUE)
    if (!length(roxygen_lines)) {
        stop(sprintf("No Roxygen2 source lines found in %s.", source_path))
    }
    if (any(grepl("^#'[[:space:]]*@param[[:space:]]*$", roxygen_lines))) {
        stop(sprintf("Empty @param tag found in %s.", source_path))
    }
    if (any(grepl("^#'[[:space:]]*@(name|rdname|aliases)[[:space:]]*$",
                  roxygen_lines))) {
        stop(sprintf("Empty topic-link tag found in %s.", source_path))
    }
    parsed <- parse(file = source_path, keep.source = FALSE)
    function_count <- sum(vapply(
        as.list(parsed),
        function(expr) {
            is.call(expr) &&
                identical(as.character(expr[[1L]]), "<-") &&
                length(expr) >= 3L &&
                is.call(expr[[3L]]) &&
                identical(as.character(expr[[3L]][[1L]]), "function")
        },
        logical(1)
    ))
    param_count <- sum(grepl("^#'[[:space:]]*@param[[:space:]]+",
                             roxygen_lines))
    export_count <- sum(grepl("^#'[[:space:]]*@export([[:space:]]|$)",
                              roxygen_lines))
    cat(sprintf(
        "SOURCE=%s ROXYGEN_LINES=%d PARAM_TAGS=%d EXPORT_TAGS=%d TOP_LEVEL_ARROW_FUNCTIONS=%d\n",
        source_path, length(roxygen_lines), param_count, export_count,
        function_count
    ))
}
RSCRIPT

git --no-optional-locks diff --check -- "${scope_paths[@]}"
git --no-optional-locks diff --cached --check -- "${scope_paths[@]}"

printf 'Documentation sources parse successfully, syntax inventory is recorded, relation output is generated-file-free, and changes since baseline are limited to Roxygen2 lines.\n'
```

## Repair Documentation Sources

```toml donna
id = "repair_documentation_sources"
kind = "donna.lib.request_action"
```

Documentation-source verification failed.

Standard output:

```text
{{ donna.lib.task_variable("source_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("source_stderr") }}
```

1. If a scoped Roxygen2 block is malformed or inconsistent, repair only `#'` lines in the authorized source and `{{ donna.lib.goto("verify_documentation_sources") }}`.
2. If the selected sources or plan are wrong and every scoped source still matches its saved baseline, update only the session plan and `{{ donna.lib.goto("validate_documentation_scope") }}` before any source edit.
3. If the workflow verifier is defective, query its relations, repair it, validate the control flow, and `{{ donna.lib.goto("verify_workflow_artifact") }}`.
4. If a scope change is needed after any scoped source differs from its saved baseline, or if HEAD, the Git index, protected worktree state, executable code, ordinary comments, generated output, or unrelated files changed, preserve the change and `{{ donna.lib.goto("record_blocker") }}`.
5. Do not weaken the Roxygen-only boundary to make verification pass.

## Review Documentation Semantics

```toml donna
id = "review_documentation_semantics"
kind = "donna.lib.request_action"
```

Deterministic documentation-source checks passed.

```text
{{ donna.lib.task_variable("source_stdout") }}
```

1. Inspect each complete scoped source and its exact Roxygen-only change since the saved baseline. If there is no change, confirm the existing documentation was still reviewed meaningfully.
2. Confirm descriptions use current terminology; every applicable formal and default is documented directly or through valid inheritance; return values, side effects, errors, deprecations, examples, and dataset fields are accurate; and shared topics are intentional.
3. Confirm export, S3, S4, class, generic, method, dataset, import, and registration directives are consistent with the editable interface source. Do not inspect generated namespace or manual output.
4. Compare mapped tests only as read-only behavioral evidence. Do not execute tests in this workflow.
5. Confirm no tool or command used by this workflow read, enumerated, compared, changed, or generated `man/` or `NAMESPACE`, and no documentation-generation, package-wide, network, staging, commit, or session-reset action occurred.
6. If a source-level repair is needed, `{{ donna.lib.goto("repair_documentation_sources") }}`.
7. If the semantic review is complete, `{{ donna.lib.goto("capture_review_snapshot") }}`.
8. If a correct conclusion requires generated output, behavior changes, or developer direction, `{{ donna.lib.goto("record_blocker") }}`.

## Capture Review Snapshot

```toml donna
id = "capture_review_snapshot"
kind = "donna.lib.run_script"
save_stdout_to = "snapshot_stdout"
save_stderr_to = "snapshot_stderr"
goto_on_success = "final_review_and_handoff"
goto_on_failure = "repair_review_snapshot"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/review-documentation-sources"
scope_file="$plan_dir/scope.txt"
baseline_pointer="$plan_dir/baseline-current"

if [[ ! -f "$scope_file" ]] || [[ -L "$scope_file" ]] ||
   [[ ! -f "$baseline_pointer" ]] || [[ -L "$baseline_pointer" ]]; then
    printf 'Scope or baseline pointer is missing or unsafe.\n' >&2
    exit 1
fi

attempt_name="$(sed -n '1p' "$baseline_pointer")"
attempt_dir="$plan_dir/$attempt_name"
if [[ ! "$attempt_name" =~ ^baseline[.][A-Za-z0-9]+$ ]] ||
   [[ ! -d "$attempt_dir" ]] || [[ -L "$attempt_dir" ]] ||
   [[ ! -f "$attempt_dir/state" ]] || [[ -L "$attempt_dir/state" ]] ||
   [[ ! -f "$attempt_dir/source-baseline.tsv" ]] ||
   [[ -L "$attempt_dir/source-baseline.tsv" ]]; then
    printf 'Current baseline attempt is missing or malformed.\n' >&2
    exit 1
fi

scope_sha="$(sha256sum < "$scope_file" | awk '{ print $1 }')"
saved_scope_sha="$(sed -n 's/^scope=//p' "$attempt_dir/state")"
if [[ "$scope_sha" != "$saved_scope_sha" ]]; then
    printf 'Scope changed before review snapshot capture.\n' >&2
    exit 1
fi

while IFS= read -r artifact_id; do
    source_path="${artifact_id#@/}"
    baseline_path="$attempt_dir/$source_path"
    if [[ ! -f "$source_path" ]] || [[ -L "$source_path" ]] ||
       [[ ! -f "$baseline_path" ]] || [[ -L "$baseline_path" ]]; then
        printf 'Scoped source or saved copy is missing or unsafe: %s\n' \
            "$source_path" >&2
        exit 1
    fi

    expected_record="$(
        awk -F '\t' -v path="$source_path" '$1 == path { print; exit }' \
            "$attempt_dir/source-baseline.tsv"
    )"
    if [[ -z "$expected_record" ]]; then
        printf 'No source baseline record exists for %s.\n' \
            "$source_path" >&2
        exit 1
    fi
    IFS=$'\t' read -r _ expected_mode expected_projection \
        <<< "$expected_record"
    actual_mode="$(stat -c '%f' -- "$source_path")"
    actual_projection="$(
        awk '!/^#'\''/' "$source_path" |
            sha256sum |
            awk '{ print $1 }'
    )"
    if [[ "$actual_mode" != "$expected_mode" ]] ||
       [[ "$actual_projection" != "$expected_projection" ]]; then
        printf 'Semantic review changed non-Roxygen content or mode in %s.\n' \
            "$source_path" >&2
        exit 1
    fi

    diff_output="$(mktemp "$plan_dir/snapshot-diff.XXXXXXXX")"
    diff_status=0
    diff -u -- "$baseline_path" "$source_path" > "$diff_output" ||
        diff_status=$?
    if (( diff_status > 1 )); then
        printf 'Could not compare source baseline for %s.\n' \
            "$source_path" >&2
        rm -f -- "$diff_output"
        exit "$diff_status"
    fi
    while IFS= read -r diff_line; do
        case "$diff_line" in
            '--- '*|'+++ '*|'@@ '*) ;;
            +*|-*)
                changed_line="${diff_line:1}"
                if [[ "$changed_line" != "#'"* ]]; then
                    printf 'Semantic review changed a non-Roxygen line in %s: %s\n' \
                        "$source_path" "$diff_line" >&2
                    rm -f -- "$diff_output"
                    exit 1
                fi
                ;;
        esac
    done < "$diff_output"
    rm -f -- "$diff_output"
done < "$scope_file"

review_fingerprint="$(
    while IFS= read -r artifact_id; do
        source_path="${artifact_id#@/}"
        printf '%s\0%s\0' "$source_path" "$(stat -c '%f' -- "$source_path")"
        sha256sum < "$source_path" |
            awk '{ printf "%s%c", $1, 0 }'
    done < "$scope_file" |
        sha256sum |
        awk '{ print $1 }'
)"
if [[ ! "$review_fingerprint" =~ ^[0-9a-f]{64}$ ]]; then
    printf 'Review snapshot did not produce a SHA-256 fingerprint.\n' >&2
    exit 1
fi

review_path="$attempt_dir/review.sha256"
if [[ -L "$review_path" ]]; then
    printf 'Review snapshot destination must not be a symlink.\n' >&2
    exit 1
fi
review_temp="$(mktemp "$attempt_dir/review.XXXXXXXX")"
printf '%s\n' "$review_fingerprint" > "$review_temp"
mv -f -- "$review_temp" "$review_path"
printf 'BASELINE_ATTEMPT=%s\nREVIEW_SHA256=%s\nSCOPE_SHA256=%s\n' \
    "$attempt_name" "$review_fingerprint" "$scope_sha"
```

## Repair Review Snapshot

```toml donna
id = "repair_review_snapshot"
kind = "donna.lib.request_action"
```

Review snapshot capture failed.

Standard error:

```text
{{ donna.lib.task_variable("snapshot_stderr") }}
```

1. If the source or scope changed intentionally after semantic review, `{{ donna.lib.goto("verify_documentation_sources") }}` and repeat the review.
2. If the snapshot implementation is defective, query the workflow relations, repair it, validate the control flow, and `{{ donna.lib.goto("verify_workflow_artifact") }}`.
3. If the state changed concurrently or the correct snapshot cannot be captured without broader authority, `{{ donna.lib.goto("record_blocker") }}`.

## Final Review and Handoff

```toml donna
id = "final_review_and_handoff"
kind = "donna.lib.request_action"
```

The reviewed source snapshot is:

```text
{{ donna.lib.task_variable("snapshot_stdout") }}
```

1. Inspect the exact scoped Roxygen-only diff or confirm the reviewed no-change outcome, the session plan, deterministic source evidence, and the complete workflow source without modifying the Git index.
2. Confirm the selected sources are ready for the maintainer to regenerate documentation and `NAMESPACE` manually with devtools. Do not perform or inspect that generation in this workflow.
3. Confirm the final report will name every reviewed source, any Roxygen changes, mapped behavioral evidence, the parse and source-integrity results, and the maintainer-owned generation handoff.
4. While Donna is awaiting this action and no workflow operation is executing, run `donna -p llm validate @/workflows/review-documentation-sources.donna.md`, `donna -p llm render @/workflows/review-documentation-sources.donna.md --mode view`, and `donna -p llm list`; inspect the rendered transitions and listed summary.
5. Requery this workflow in both governance directions and confirm the catalog, agent guidance, scalable Depmesh rules, package exclusion, and artifact form one coherent implementation.
6. If the review and handoff are complete, `{{ donna.lib.goto("verify_scope_stability") }}` so protected state is checked immediately before finish.
7. If a source-level discrepancy remains, `{{ donna.lib.goto("repair_documentation_sources") }}`.
8. If a contract conflict, generated-output dependency, or ownership decision remains, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Scope Stability

```toml donna
id = "verify_scope_stability"
kind = "donna.lib.run_script"
save_stdout_to = "stability_stdout"
save_stderr_to = "stability_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_scope_stability"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

plan_dir=".session/donna/review-documentation-sources"
scope_file="$plan_dir/scope.txt"
baseline_pointer="$plan_dir/baseline-current"

if [[ ! -f "$scope_file" ]] || [[ -L "$scope_file" ]] ||
   [[ ! -f "$baseline_pointer" ]] || [[ -L "$baseline_pointer" ]]; then
    printf 'Scope or baseline pointer is missing or unsafe.\n' >&2
    exit 1
fi

attempt_name="$(sed -n '1p' "$baseline_pointer")"
attempt_dir="$plan_dir/$attempt_name"
if [[ ! "$attempt_name" =~ ^baseline[.][A-Za-z0-9]+$ ]] ||
   [[ ! -d "$attempt_dir" ]] || [[ -L "$attempt_dir" ]] ||
   [[ ! -f "$attempt_dir/state" ]] || [[ -L "$attempt_dir/state" ]] ||
   [[ ! -f "$attempt_dir/review.sha256" ]] ||
   [[ -L "$attempt_dir/review.sha256" ]]; then
    printf 'Current baseline or review snapshot is missing or malformed.\n' >&2
    exit 1
fi

state_value() {
    local key="$1"
    sed -n "s/^${key}=//p" "$attempt_dir/state"
}

emit_worktree_manifest() {
    git --no-optional-locks ls-files -z --cached --others --exclude-standard -- "$@" |
        LC_ALL=C sort -z |
        while IFS= read -r -d '' path; do
            printf '%s\0' "$path" || exit $?
            if [[ -L "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'symlink\0%s\0' "$mode" || exit $?
                readlink -z -- "$path" |
                    sha256sum |
                    awk '{ printf "%s%c", $1, 0 }'
            elif [[ -f "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'file\0%s\0' "$mode" || exit $?
                sha256sum < "$path" |
                    awk '{ printf "%s%c", $1, 0 }'
            elif [[ -e "$path" ]]; then
                mode="$(stat -c '%f' -- "$path")"
                printf 'other\0%s\0-\0' "$mode" || exit $?
            else
                printf 'missing\0-\0-\0' || exit $?
            fi
        done
}

emit_index_manifest() {
    printf 'index-stage-v1\0'
    git --no-optional-locks ls-files --stage -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-assume-skip-v1\0'
    git --no-optional-locks ls-files --stage -v -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-fsmonitor-v1\0'
    git --no-optional-locks ls-files --stage -f -z -- "$@" |
        LC_ALL=C sort -z
    printf '\0index-entry-flags-v1\0'
    LC_ALL=C git --no-optional-locks ls-files --stage --debug -z -- "$@"
    printf '\0index-resolve-undo-v1\0'
    git --no-optional-locks ls-files --resolve-undo -z -- "$@" |
        LC_ALL=C sort -z
}

scope_sha="$(sha256sum < "$scope_file" | awk '{ print $1 }')"
if [[ "$scope_sha" != "$(state_value scope)" ]]; then
    printf 'Scope changed after baseline capture.\n' >&2
    exit 1
fi

scope_paths=()
while IFS= read -r artifact_id; do
    scope_paths+=("${artifact_id#@/}")
done < "$scope_file"

protected_pathspecs=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
for source_path in "${scope_paths[@]}"; do
    protected_pathspecs+=(":(top,exclude)${source_path}")
done

current_head="$(git --no-optional-locks rev-parse --verify 'HEAD^{commit}')"
current_index="$(
    emit_index_manifest \
        . \
        ':(top,exclude)man/**' \
        ':(top,exclude)NAMESPACE' |
        sha256sum |
        awk '{ print $1 }'
)"
current_protected="$(
    {
        git --no-optional-locks status --porcelain=v2 -z \
            --untracked-files=all -- "${protected_pathspecs[@]}" &&
            printf '\0protected-worktree-v1\0' &&
            emit_worktree_manifest "${protected_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"
current_review="$(
    while IFS= read -r artifact_id; do
        source_path="${artifact_id#@/}"
        if [[ ! -f "$source_path" ]] || [[ -L "$source_path" ]]; then
            printf 'Reviewed source is missing or unsafe: %s\n' \
                "$source_path" >&2
            exit 1
        fi
        printf '%s\0%s\0' "$source_path" "$(stat -c '%f' -- "$source_path")"
        sha256sum < "$source_path" |
            awk '{ printf "%s%c", $1, 0 }'
    done < "$scope_file" |
        sha256sum |
        awk '{ print $1 }'
)"
saved_review="$(sed -n '1p' "$attempt_dir/review.sha256")"

if [[ "$current_head" != "$(state_value head)" ]]; then
    printf 'HEAD changed during documentation-source review.\n' >&2
    exit 1
fi
if [[ "$current_index" != "$(state_value index)" ]]; then
    printf 'The non-generated Git index changed during documentation-source review.\n' >&2
    exit 1
fi
if [[ "$current_protected" != "$(state_value protected)" ]]; then
    printf 'A protected non-generated worktree artifact changed during documentation-source review.\n' >&2
    exit 1
fi
if [[ "$current_review" != "$saved_review" ]]; then
    printf 'A reviewed documentation source changed after snapshot capture.\n' >&2
    exit 1
fi

safe_diff_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks diff --check -- "${safe_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_diff_paths[@]}"

workflow_id="@/workflows/review-documentation-sources.donna.md"
governance_output="$(depmesh -p llm dependencies "$workflow_id")"
for governing_spec in \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
do
    if ! printf '%s\n' "$governance_output" |
        rg -Fqx -- "- ${governing_spec}"; then
        printf 'Final governance query is missing %s.\n' "$governing_spec" >&2
        exit 1
    fi
    if ! depmesh -p llm dependencies --relation governs "$governing_spec" |
        sed -n 's/^- \(@\/.*\)$/\1/p' |
        rg -Fqx -- "$workflow_id"; then
        printf 'Final reverse governance query is missing %s for %s.\n' \
            "$workflow_id" "$governing_spec" >&2
        exit 1
    fi
done

printf 'HEAD=%s\nINDEX_SHA256=%s\nPROTECTED_SHA256=%s\nREVIEW_SHA256=%s\nSCOPE_SHA256=%s\n' \
    "$current_head" "$current_index" "$current_protected" \
    "$current_review" "$scope_sha"
printf 'Reviewed documentation sources and protected non-generated project state match their approved snapshots.\n'
```

## Repair Scope Stability

```toml donna
id = "repair_scope_stability"
kind = "donna.lib.request_action"
```

The final protected-scope comparison failed.

Standard output:

```text
{{ donna.lib.task_variable("stability_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("stability_stderr") }}
```

1. Inspect only explicitly allowed non-generated status, scoped source diffs, and session evidence without changing HEAD or the Git index.
2. If the reviewed source changed intentionally, `{{ donna.lib.goto("verify_documentation_sources") }}` and repeat semantic review and snapshot capture.
3. If protected state changed concurrently or outside the authorized scope, preserve it and `{{ donna.lib.goto("record_blocker") }}`.
4. If the comparator is defective, query the workflow relations, repair only the verifier, validate it, and `{{ donna.lib.goto("verify_workflow_artifact") }}`.
5. Do not inspect generated output to diagnose or bypass this gate.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The selected public API documentation sources passed scoped Depmesh discovery, R parsing, Roxygen-only edit enforcement, semantic review, workflow integration checks, and protected-state verification. Report the reviewed sources, Roxygen changes or deliberate no-change result, mapped behavioral evidence, validation and rendering evidence, and the explicit handoff for the maintainer to regenerate documentation and `NAMESPACE` manually. Generated output remains uninspected.

If this was the first representative execution and the catalog still says `implementation in progress`, mark it `implemented` only after this finish. Then rerun this workflow in the final catalog state and repeat workflow validation, view rendering, listing, bidirectional governance queries, and non-generated diff checks before reporting completion.

## Record Blocker

```toml donna
id = "record_blocker"
kind = "donna.lib.request_action"
```

1. Preserve the exact completed gates, scoped sources, failure evidence, and external action or decision required.
2. Determine whether another cataloged workflow owns the resolution. If so, name that workflow and stage without running it unless separately authorized or required.
3. If no remaining workflow resolves the finding, query all relations for `{{ donna.lib.path("@/workflows/reports/documentation-source-open-questions.md") }}` and create or update that exact report with the source path, public interface, evidence, ownership question, and required decision. Do not use the `.donna.md` extension or add the report to the workflow catalog.
4. Preserve unrelated work and do not inspect generated output, change behavior, stage files, create commits, use the network, or reset the Donna session.
5. After recording the blocker, `{{ donna.lib.goto("blocked") }}`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The documentation-source workflow cannot complete without developer or maintainer action. Report completed gates, reviewed sources, the exact source-level, generated-output, environment, contract, or authority blocker, any persistent open-question entry, and the decision required. Do not inspect generated output, change package behavior, stage files, create commits, use the network, or reset the Donna session.
