# Verify Specifications

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "capture_scope_baseline"
```

Verify structural and semantic consistency across every project specification, its index, and its Depmesh governance without changing package implementation behavior.

This workflow MUST NOT access `man/`, use the network, change package behavior, stage files, create commits, or reset the Donna session.
Every Git inspection MUST use `git --no-optional-locks` so read-only checks do not persist optional index refreshes.

## Capture Scope Baseline

```toml donna
id = "capture_scope_baseline"
kind = "donna.lib.run_script"
save_stdout_to = "scope_baseline"
save_stderr_to = "scope_baseline_stderr"
goto_on_success = "preflight"
goto_on_failure = "repair_scope_baseline"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

emit_worktree_manifest() {
    git --no-optional-locks ls-files -z --cached --others --exclude-standard -- "$@" |
        LC_ALL=C sort -z |
        while IFS= read -r -d '' path; do
            printf '%s\0' "$path" || exit $?
            if [[ -L "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect symlink mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'symlink\0%s\0' "$mode" || exit $?
                if readlink -z -- "$path" |
                    sha256sum |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint symlink target: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -f "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect file mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'file\0%s\0' "$mode" || exit $?
                if sha256sum < "$path" |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint file contents: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -e "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect artifact mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'other\0%s\0-\0' "$mode" || exit $?
            else
                printf 'missing\0-\0-\0' || exit $?
            fi
        done
}

emit_index_manifest() {
    printf 'index-stage-v3\0' || return
    git --no-optional-locks ls-files --stage -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-assume-skip-v1\0' || return
    git --no-optional-locks ls-files --stage -v -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-fsmonitor-v1\0' || return
    git --no-optional-locks ls-files --stage -f -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-entry-flags-v1\0' || return
    LC_ALL=C git --no-optional-locks ls-files --stage --debug -z -- "$@" || return
    printf '\0index-resolve-undo-v1\0' || return
    git --no-optional-locks ls-files --resolve-undo -z -- "$@" |
        LC_ALL=C sort -z
}

existing_baseline="$(
    cat <<'DONNA_SCOPE_BASELINE'
{{ donna.lib.task_variable("scope_baseline") }}
DONNA_SCOPE_BASELINE
)"
if [[ "$existing_baseline" =~ ^v3:[0-9a-f]{40,64}:[0-9a-f]{64}:[0-9a-f]{64}:[0-9a-f]{64}$ ]]; then
    printf '%s\n' "$existing_baseline"
    exit 0
fi
if [[ -n "$existing_baseline" ]] &&
   { [[ "$existing_baseline" != '$$donna task_variable variable '* ]] ||
     [[ "$existing_baseline" != *"does not found."* ]]; }; then
    printf 'Existing scope baseline has an invalid format; refusing to replace it.\n' >&2
    exit 1
fi

protected_pathspecs=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)AGENTS.md'
    ':(top,exclude)depmesh.toml'
    ':(top,exclude)specs/**'
    ':(top,exclude)workflows/verify-specifications.donna.md'
)

index_fingerprint="$(
    emit_index_manifest \
        . \
        ':(top,exclude)man/**' |
        sha256sum |
        awk '{ print $1 }'
)"

protected_worktree_fingerprint="$(
    {
        git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
            "${protected_pathspecs[@]}" &&
            printf '\0protected-worktree-manifest-v1\0' &&
            emit_worktree_manifest "${protected_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"

review_pathspecs=(
    AGENTS.md
    depmesh.toml
    specs
    workflows/verify-specifications.donna.md
)
initial_review_fingerprint="$(
    {
        git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
            "${review_pathspecs[@]}" &&
            printf '\0review-worktree-manifest-v1\0' &&
            emit_worktree_manifest "${review_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"

head_revision="$(git --no-optional-locks rev-parse --verify 'HEAD^{commit}')"
if [[ ! "$head_revision" =~ ^[0-9a-f]{40,64}$ ]] ||
   [[ ! "$index_fingerprint" =~ ^[0-9a-f]{64}$ ]] ||
   [[ ! "$protected_worktree_fingerprint" =~ ^[0-9a-f]{64}$ ]] ||
   [[ ! "$initial_review_fingerprint" =~ ^[0-9a-f]{64}$ ]]; then
    printf 'Scope baseline did not produce the expected Git and SHA-256 identifiers.\n' >&2
    exit 1
fi

printf 'v3:%s:%s:%s:%s\n' \
    "$head_revision" \
    "$index_fingerprint" \
    "$protected_worktree_fingerprint" \
    "$initial_review_fingerprint"
```

## Repair Scope Baseline

```toml donna
id = "repair_scope_baseline"
kind = "donna.lib.request_action"
```

The workflow could not capture its initial protected-scope baseline.

Standard error:

```text
{{ donna.lib.task_variable("scope_baseline_stderr") }}
```

1. Diagnose only the baseline command, repository state, required local command availability, or an incorrect verifier expectation.
2. Do not change HEAD, the Git index, package implementation, unrelated files, or any path under `man/`.
3. If a verifier repair is sufficient, `{{ donna.lib.goto("capture_scope_baseline") }}`.
4. If an initial non-`man/` baseline cannot be captured without developer direction, `{{ donna.lib.goto("blocked") }}`.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/meta/general.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}` completely.
2. Read every other Markdown specification under `{{ donna.lib.path("@/specs") }}` completely. Do not infer specification content from the index summary.
3. Read the built-in Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations separately for every specification, `{{ donna.lib.path("@/AGENTS.md") }}`, this workflow, and each other artifact already known to require an edit.
4. Inspect the current worktree with `git --no-optional-locks` without modifying the Git index. Preserve unrelated user changes.
5. Confirm that this action belongs to the current Donna task and that no unrelated action request is being displaced. Do not reset the session.
6. Confirm that this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna execution process was active.
7. Do not read, enumerate, compare, or otherwise inspect any path under `man/`.
8. If the specification hierarchy, current diff, and applicable relations are understood, `{{ donna.lib.goto("run_index_and_structure") }}`.
9. If required context is unavailable or the specifications already conflict about how this workflow must operate, `{{ donna.lib.goto("blocked") }}`.

## Restart Specification Verification

```toml donna
id = "verify_index_and_structure"
kind = "donna.lib.output"
next_operation_id = "capture_scope_baseline"
```

Preserve the original protected-scope baseline, repeat preflight after the repair, and restart every deterministic and semantic specification gate.

## Run Index and Structure

```toml donna
id = "run_index_and_structure"
kind = "donna.lib.run_script"
save_stdout_to = "structure_stdout"
save_stderr_to = "structure_stderr"
goto_on_success = "verify_definite_path_style"
goto_on_failure = "repair_index_and_structure"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

count_arguments() {
    printf '%d\n' "$#"
}

actual_specs_output="$(
    find specs -type f -name '*.md' -print |
        LC_ALL=C sort
)"
actual_specs=()
if [[ -n "$actual_specs_output" ]]; then
    mapfile -t actual_specs <<< "$actual_specs_output"
fi
actual_spec_count="$(count_arguments "${actual_specs[@]}")"
if (( actual_spec_count == 0 )); then
    printf 'No Markdown specifications were found under specs/.\n' >&2
    exit 1
fi

for required_index_section in \
    "Specification directories" \
    "Specification documents"
do
    section_count="$(
        awk -v heading="## ${required_index_section}" '
            $0 == heading {
                count++
            }
            END {
                print count + 0
            }
        ' specs/intro.md
    )"
    if (( section_count != 1 )); then
        printf 'Specification index must contain exactly one %s section; found %d.\n' \
            "$required_index_section" "$section_count" >&2
        exit 1
    fi
done

indexed_specs_output="$(
    awk '
        $0 == "## Specification documents" {
            in_section = 1
            next
        }
        in_section && /^## / {
            exit
        }
        in_section && /^- `/ {
            line = $0
            if (line !~ /^- `\.\/specs\/[^`]+\.md`[[:space:]]+[^[:space:]]/) {
                next
            }
            sub(/^- `/, "", line)
            sub(/`.*/, "", line)
            print line
        }
    ' specs/intro.md |
        LC_ALL=C sort
)"
indexed_specs=()
if [[ -n "$indexed_specs_output" ]]; then
    mapfile -t indexed_specs <<< "$indexed_specs_output"
fi
document_bullet_count="$(
    awk '
        $0 == "## Specification documents" {
            in_section = 1
            next
        }
        in_section && /^## / {
            exit
        }
        in_section && /^- / {
            count++
        }
        END {
            print count + 0
        }
    ' specs/intro.md
)"
indexed_spec_count="$(count_arguments "${indexed_specs[@]}")"
if (( indexed_spec_count != document_bullet_count )); then
    printf 'Every specification-document bullet must begin with a `./specs/...md` path.\n' >&2
    exit 1
fi

duplicate_specs="$(
    printf '%s\n' "${indexed_specs[@]}" |
        uniq -d
)"
if [[ -n "$duplicate_specs" ]]; then
    printf 'Specification index contains duplicate document entries:\n%s\n' \
        "$duplicate_specs" >&2
    exit 1
fi

actual_spec_entries="$(
    printf './%s\n' "${actual_specs[@]}" |
        LC_ALL=C sort
)"
indexed_spec_entries="$(
    printf '%s\n' "${indexed_specs[@]}" |
        LC_ALL=C sort
)"
if [[ "$actual_spec_entries" != "$indexed_spec_entries" ]]; then
    printf 'Specification files and index entries differ.\nActual:\n%s\nIndexed:\n%s\n' \
        "$actual_spec_entries" "$indexed_spec_entries" >&2
    exit 1
fi

actual_directories_output="$(
    find specs -type d -print |
        while IFS= read -r directory; do
            if [[ "$directory" == "specs" ]]; then
                printf './specs/\n'
            else
                printf './%s/\n' "$directory"
            fi
        done |
        LC_ALL=C sort
)"
actual_directories=()
if [[ -n "$actual_directories_output" ]]; then
    mapfile -t actual_directories <<< "$actual_directories_output"
fi
indexed_directories_output="$(
    awk '
        $0 == "## Specification directories" {
            in_section = 1
            next
        }
        in_section && /^## / {
            exit
        }
        in_section && /^- `/ {
            line = $0
            if (line !~ /^- `\.\/specs[^`]*\/`[[:space:]]+[^[:space:]]/) {
                next
            }
            sub(/^- `/, "", line)
            sub(/`.*/, "", line)
            print line
        }
    ' specs/intro.md |
        LC_ALL=C sort
)"
indexed_directories=()
if [[ -n "$indexed_directories_output" ]]; then
    mapfile -t indexed_directories <<< "$indexed_directories_output"
fi
directory_bullet_count="$(
    awk '
        $0 == "## Specification directories" {
            in_section = 1
            next
        }
        in_section && /^## / {
            exit
        }
        in_section && /^- / {
            count++
        }
        END {
            print count + 0
        }
    ' specs/intro.md
)"
indexed_directory_count="$(count_arguments "${indexed_directories[@]}")"
if (( indexed_directory_count != directory_bullet_count )); then
    printf 'Every specification-directory bullet must contain a `./specs.../` path.\n' >&2
    exit 1
fi

duplicate_directories="$(
    printf '%s\n' "${indexed_directories[@]}" |
        uniq -d
)"
if [[ -n "$duplicate_directories" ]]; then
    printf 'Specification index contains duplicate directory entries:\n%s\n' \
        "$duplicate_directories" >&2
    exit 1
fi

actual_directory_entries="$(
    printf '%s\n' "${actual_directories[@]}" |
        LC_ALL=C sort
)"
indexed_directory_entries="$(
    printf '%s\n' "${indexed_directories[@]}" |
        LC_ALL=C sort
)"
if [[ "$actual_directory_entries" != "$indexed_directory_entries" ]]; then
    printf 'Specification directories and index entries differ.\nActual:\n%s\nIndexed:\n%s\n' \
        "$actual_directory_entries" "$indexed_directory_entries" >&2
    exit 1
fi

markdown_headings() {
    local spec_path="$1"
    awk '
        {
            normalized = $0
            indent = 0
            while (indent < 3 && substr(normalized, 1, 1) == " ") {
                normalized = substr(normalized, 2)
                indent++
            }

            fence_character = substr(normalized, 1, 1)
            fence_length = 0
            if (fence_character == "`" || fence_character == "~") {
                while (substr(normalized, fence_length + 1, 1) == fence_character) {
                    fence_length++
                }
            }
            fence_remainder = substr(normalized, fence_length + 1)

            if (in_fence) {
                if (fence_character == active_fence_character &&
                    fence_length >= active_fence_length &&
                    fence_remainder ~ /^[[:space:]]*$/) {
                    in_fence = 0
                    active_fence_character = ""
                    active_fence_length = 0
                }
                next
            }
            if (fence_length >= 3) {
                in_fence = 1
                active_fence_character = fence_character
                active_fence_length = fence_length
                previous_text = ""
                previous_line = 0
                next
            }

            if (normalized ~ /^# / || normalized ~ /^## /) {
                heading = normalized
                sub(/[[:space:]]+#+[[:space:]]*$/, "", heading)
                print NR "\t" heading
                previous_text = ""
                previous_line = 0
                next
            }

            if (normalized ~ /^=+[[:space:]]*$/ && previous_text != "") {
                print previous_line "\t# " previous_text
                previous_text = ""
                previous_line = 0
                next
            }
            if (normalized ~ /^-+[[:space:]]*$/ && previous_text != "") {
                print previous_line "\t## " previous_text
                previous_text = ""
                previous_line = 0
                next
            }

            candidate = normalized
            sub(/[[:space:]]+$/, "", candidate)
            if (candidate != "" &&
                candidate !~ /^([#>]|[-*+][[:space:]]|[0-9]+[.)][[:space:]])/) {
                previous_text = candidate
                previous_line = NR
            } else {
                previous_text = ""
                previous_line = 0
            }
        }
    ' "$spec_path"
}

for spec_path in "${actual_specs[@]}"; do
    heading_records_output="$(markdown_headings "$spec_path")"
    heading_records=()
    if [[ -n "$heading_records_output" ]]; then
        mapfile -t heading_records <<< "$heading_records_output"
    fi
    h1_count=0
    h1_line=""
    h1_title=""
    h2_lines=()
    h2_titles=()
    for heading_record in "${heading_records[@]}"; do
        line_number="${heading_record%%$'\t'*}"
        heading="${heading_record#*$'\t'}"
        if [[ "$heading" == "# "* ]]; then
            h1_count=$((h1_count + 1))
            h1_line="$line_number"
            h1_title="${heading#\# }"
        elif [[ "$heading" == "## "* ]]; then
            h2_lines+=("$line_number")
            h2_titles+=("${heading#\#\# }")
        fi
    done

    if (( h1_count != 1 )); then
        printf '%s must contain exactly one H1 heading; found %d.\n' \
            "$spec_path" "$h1_count" >&2
        exit 1
    fi

    first_h2_line="${h2_lines[0]:-}"
    if [[ -z "$first_h2_line" ]] || (( h1_line >= first_h2_line )); then
        printf '%s must place its single H1 before its top-level sections.\n' \
            "$spec_path" >&2
        exit 1
    fi

    if [[ -z "$h1_title" ]]; then
        printf '%s must give its H1 heading a nonempty title.\n' \
            "$spec_path" >&2
        exit 1
    fi

    h2_count="$(count_arguments "${h2_titles[@]}")"
    if (( h2_count < 2 )); then
        printf '%s must contain Goal of the document and Scope sections.\n' \
            "$spec_path" >&2
        exit 1
    fi
    if [[ "${h2_titles[0]}" != "Goal of the document" ]] ||
       [[ "${h2_titles[1]}" != "Scope" ]]; then
        printf '%s must begin with Goal of the document followed by Scope.\n' \
            "$spec_path" >&2
        exit 1
    fi

    dictionary_count=0
    for h2_title in "${h2_titles[@]}"; do
        if [[ "$h2_title" == "Dictionary" ]]; then
            dictionary_count=$((dictionary_count + 1))
        fi
    done
    if (( dictionary_count > 1 )); then
        printf '%s contains more than one Dictionary section.\n' \
            "$spec_path" >&2
        exit 1
    fi
    if (( dictionary_count == 1 )) &&
       { (( h2_count < 3 )) || [[ "${h2_titles[2]}" != "Dictionary" ]]; }; then
        printf '%s must place Dictionary immediately after Scope.\n' \
            "$spec_path" >&2
        exit 1
    fi
done

actual_directory_count="$(count_arguments "${actual_directories[@]}")"
printf 'Matched %d specification files and %d directories to unique index entries; every specification has one H1 and the required opening H2 order outside fenced code blocks.\n' \
    "$actual_spec_count" "$actual_directory_count"
```

## Repair Index and Structure

```toml donna
id = "repair_index_and_structure"
kind = "donna.lib.request_action"
```

The specification index or document-structure gate failed.

Standard output:

```text
{{ donna.lib.task_variable("structure_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("structure_stderr") }}
```

1. Query all configured relations separately for every artifact to be edited.
2. Repair only the specification index, specification structure, affected scalable governance rules, affected agent guidance, or an incorrect verifier expectation.
3. When a specification is added, moved, renamed, or removed, synchronize `{{ donna.lib.path("@/specs/intro.md") }}` and every affected governance reference in the same repair.
4. Preserve the intended contract, unrelated user changes, and the Git index. Do not change package implementation behavior or inspect `man/`.
5. After repair, `{{ donna.lib.goto("verify_index_and_structure") }}`.
6. If the required structure conflicts with another specification or a path/ownership decision is needed, `{{ donna.lib.goto("blocked") }}`.

## Verify Definite Path Style

```toml donna
id = "verify_definite_path_style"
kind = "donna.lib.run_script"
save_stdout_to = "path_style_stdout"
save_stderr_to = "path_style_stderr"
goto_on_success = "review_specification_semantics"
goto_on_failure = "repair_path_style"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

spec_paths_output="$(
    find specs -type f -name '*.md' -print |
        LC_ALL=C sort
)"
spec_paths=()
if [[ -n "$spec_paths_output" ]]; then
    mapfile -t spec_paths <<< "$spec_paths_output"
fi
if [[ -z "$spec_paths_output" ]]; then
    printf 'No Markdown specifications were found for path-style verification.\n' >&2
    exit 1
fi

bare_project_path_pattern='`(?:AGENTS\.md|DESCRIPTION|NAMESPACE|NEWS(?:\.md)?|\.Rbuildignore|depmesh\.toml|donna\.toml|(?:R|tests|vignettes|data|inst|specs|workflows|bin|man)/[^`]*)`'
if violations="$(
    rg -n --pcre2 "$bare_project_path_pattern" "${spec_paths[@]}"
)"; then
    printf 'Project-root-relative paths in specification prose must use the ./ form:\n%s\n' \
        "$violations" >&2
    exit 1
else
    status=$?
    if (( status != 1 )); then
        printf 'Path-style scan failed with exit code %d.\n' "$status" >&2
        exit "$status"
    fi
fi

malformed_root_id_pattern='@/(?:/|\.\.(?:/|`|[[:space:]]|$)|[^`[:space:]]+/\.\.(?:/|`|[[:space:]]|$))'
if violations="$(
    rg -n --pcre2 "$malformed_root_id_pattern" "${spec_paths[@]}"
)"; then
    printf 'Malformed root-anchored Depmesh artifact identifier:\n%s\n' \
        "$violations" >&2
    exit 1
else
    status=$?
    if (( status != 1 )); then
        printf 'Artifact-identifier scan failed with exit code %d.\n' "$status" >&2
        exit "$status"
    fi
fi

printf 'No definite bare project path or malformed root-anchored artifact identifier was found; ambiguous prose remains subject to semantic review.\n'
```

## Repair Path Style

```toml donna
id = "repair_path_style"
kind = "donna.lib.request_action"
```

The definite path-style gate failed.

Standard output:

```text
{{ donna.lib.task_variable("path_style_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("path_style_stderr") }}
```

Query all configured relations for each affected specification, then repair project-root-relative prose paths to use `./path/to/artifact` and Depmesh query inputs or outputs to use `@/path/to/artifact`. Preserve code symbols, placeholders, URLs, and command names that are not project paths.

After repair, `{{ donna.lib.goto("verify_index_and_structure") }}`. If a token's meaning requires a developer decision, `{{ donna.lib.goto("blocked") }}`.

## Review Specification Semantics

```toml donna
id = "review_specification_semantics"
kind = "donna.lib.request_action"
```

The deterministic index, structure, and definite path-style gates passed.

1. Read every specification completely, including unchanged files. Do not read generated manuals.
2. Confirm each `Goal of the document` describes the specification's content and purpose without defining requirements for the document itself.
3. Confirm each `Scope` defines boundaries, identifies important exclusions when needed, and does not redirect requirements to another specification.
4. Confirm H1 names are unique and each `Dictionary` contains only specification-specific terms and appears only when those terms are needed.
5. Inventory uppercase RFC 2119 terms and lowercase requirement-like wording. Confirm every normative statement uses `MUST`, `MUST NOT`, `SHOULD`, `SHOULD NOT`, or `MAY`; is testable or reviewable; and does not depend on an undefined vague qualifier.
6. Review every project path in prose and every Depmesh query input or output, including cases the definite-token scan cannot classify.
7. Confirm each specification owns one cohesive concern; overlapping requirements have one owner; summaries do not duplicate normative requirements; and no current normative statements conflict.
8. Confirm Markdown paragraphs are not hard-wrapped merely to fit a fixed width, requirements remain at the highest precise abstraction level, and concrete artifacts are named only when they are stable contracts.
9. Compare any changed requirement with the artifacts returned by its reverse governance query. Do not change package implementation behavior in this workflow and do not silently normalize a specification to historical behavior.
10. If the specifications are clear, non-conflicting, and internally consistent, `{{ donna.lib.goto("verify_governance") }}`.
11. If a discrepancy can be repaired without changing intended implementation behavior, `{{ donna.lib.goto("repair_semantic_findings") }}`.
12. If the intended contract is ambiguous or conflicts with implementation and needs a developer decision, record the unresolved finding in an open-questions report under `{{ donna.lib.path("@/workflows/reports") }}` when no remaining workflow can resolve it, then `{{ donna.lib.goto("blocked") }}`.

## Repair Semantic Findings

```toml donna
id = "repair_semantic_findings"
kind = "donna.lib.request_action"
```

1. Query all configured relations separately for every specification or governed artifact to be edited.
2. Repair only wording, ownership, index summaries, path references, governance metadata, or verifier instructions that can be corrected without changing intended package behavior.
3. Keep normative requirements testable, use RFC 2119 keywords, and avoid duplicating requirements across specifications.
4. Preserve unrelated changes and leave the Git index untouched.
5. After repair, reread every changed specification and `{{ donna.lib.goto("verify_index_and_structure") }}`.
6. If a repair would choose a new contract or require package implementation changes, record an unresolved finding when required and `{{ donna.lib.goto("blocked") }}`.

## Verify Governance

```toml donna
id = "verify_governance"
kind = "donna.lib.run_script"
save_stdout_to = "governance_stdout"
save_stderr_to = "governance_stderr"
goto_on_success = "verify_workflow_artifact"
goto_on_failure = "repair_governance"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

count_arguments() {
    printf '%d\n' "$#"
}

relation_artifacts() {
    local relation="$1"
    local artifact="$2"
    local query_output malformed_bullets artifacts duplicates status
    if query_output="$(
        depmesh -p llm dependencies --relation "$relation" "$artifact"
    )"; then
        :
    else
        status=$?
        printf 'Depmesh %s query failed for %s with exit code %d.\n' \
            "$relation" "$artifact" "$status" >&2
        return "$status"
    fi
    if rg -Fqx -- '## warnings' <<< "$query_output" >/dev/null; then
        printf 'Depmesh returned warnings for %s on %s:\n%s\n' \
            "$relation" "$artifact" "$query_output" >&2
        return 1
    else
        status=$?
        if (( status != 1 )); then
            printf 'Could not scan Depmesh output for warnings on %s for %s.\n' \
                "$relation" "$artifact" >&2
            return "$status"
        fi
    fi
    if malformed_bullets="$(
        printf '%s\n' "$query_output" |
            sed -n '/^- /p' |
            sed -n '\#^- @/[^[:space:]]\+$#!p'
    )"; then
        :
    else
        status=$?
        printf 'Could not validate Depmesh artifact bullets for %s on %s.\n' \
            "$relation" "$artifact" >&2
        return "$status"
    fi
    if [[ -n "$malformed_bullets" ]]; then
        printf 'Depmesh returned malformed artifact bullets for %s on %s:\n%s\n' \
            "$relation" "$artifact" "$malformed_bullets" >&2
        return 1
    fi
    if artifacts="$(
        printf '%s\n' "$query_output" |
            sed -n 's/^- \(@\/[^[:space:]]\+\)$/\1/p' |
            LC_ALL=C sort
    )"; then
        :
    else
        status=$?
        printf 'Could not parse Depmesh artifacts for %s on %s.\n' \
            "$relation" "$artifact" >&2
        return "$status"
    fi
    if duplicates="$(
        printf '%s\n' "$artifacts" |
            sed '/^$/d' |
            uniq -d
    )"; then
        :
    else
        status=$?
        printf 'Could not scan Depmesh artifacts for duplicates for %s on %s.\n' \
            "$relation" "$artifact" >&2
        return "$status"
    fi
    if [[ -n "$duplicates" ]]; then
        printf 'Depmesh returned duplicate artifacts for %s on %s:\n%s\n' \
            "$relation" "$artifact" "$duplicates" >&2
        return 1
    fi
    printf '%s\n' "$artifacts"
}

assert_exact() {
    local relation="$1"
    local artifact="$2"
    shift 2
    local actual expected status
    if actual="$(relation_artifacts "$relation" "$artifact")"; then
        :
    else
        status=$?
        return "$status"
    fi
    if (( $# )); then
        expected="$(printf '%s\n' "$@" | LC_ALL=C sort)"
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
    local actual status
    if actual="$(relation_artifacts "$relation" "$artifact")"; then
        :
    else
        status=$?
        return "$status"
    fi
    if rg -Fqx -- "$expected" <<< "$actual" >/dev/null; then
        :
    else
        status=$?
        printf '%s result for %s does not contain %s.\n' \
            "$relation" "$artifact" "$expected" >&2
        if (( status == 1 )); then
            return 1
        fi
        return "$status"
    fi
}

project_root="$(pwd -P)"
assert_existing_project_artifact() {
    local artifact="$1"
    local relative resolved
    if [[ "$artifact" != @/* ]]; then
        printf 'Relation returned a non-root-anchored artifact: %s\n' \
            "$artifact" >&2
        return 1
    fi
    if [[ "$artifact" == "@/man" || "$artifact" == "@/man/"* ]]; then
        printf 'Relation returned a forbidden generated-manual artifact: %s\n' \
            "$artifact" >&2
        return 1
    fi
    if [[ "$artifact" == *"/../"* || "$artifact" == *"/.." ]]; then
        printf 'Relation returned parent traversal: %s\n' "$artifact" >&2
        return 1
    fi
    relative="${artifact#@/}"
    if [[ ! -f "$relative" ]]; then
        printf 'Relation returned a missing or non-file artifact: %s\n' \
            "$artifact" >&2
        return 1
    fi
    resolved="$(realpath -e -- "$relative")"
    case "$resolved" in
        "$project_root"/*) ;;
        *)
            printf 'Relation artifact resolves outside the project: %s\n' \
                "$artifact" >&2
            return 1
            ;;
    esac
}

if relation_output="$(depmesh -p llm relations)"; then
    :
else
    status=$?
    printf 'Depmesh relation discovery failed with exit code %d.\n' "$status" >&2
    exit "$status"
fi
for relation_id in governed_by governs; do
    if rg -Fqx -- "## ${relation_id}" <<< "$relation_output" >/dev/null; then
        :
    else
        status=$?
        printf 'Required governance relation is missing: %s\n' "$relation_id" >&2
        if (( status == 1 )); then
            exit 1
        fi
        exit "$status"
    fi
done

spec_paths_output="$(
    find specs -type f -name '*.md' -print |
        LC_ALL=C sort
)"
spec_paths=()
if [[ -n "$spec_paths_output" ]]; then
    mapfile -t spec_paths <<< "$spec_paths_output"
fi
if [[ -z "$spec_paths_output" ]]; then
    printf 'No Markdown specifications were found for governance verification.\n' >&2
    exit 1
fi
spec_ids=()
meta_governed_ids=("@/AGENTS.md")
for spec_path in "${spec_paths[@]}"; do
    spec_id="@/${spec_path}"
    spec_ids+=("$spec_id")

    if [[ "$spec_id" == "@/specs/meta/general.md" ]]; then
        assert_exact governed_by "$spec_id"
    else
        assert_contains governed_by "$spec_id" @/specs/meta/general.md
        meta_governed_ids+=("$spec_id")
    fi
    assert_contains governs "$spec_id" @/AGENTS.md
done

assert_exact governed_by @/AGENTS.md "${spec_ids[@]}"
assert_exact governs @/specs/meta/general.md "${meta_governed_ids[@]}"

for spec_id in "${spec_ids[@]}"; do
    governed_artifacts="$(relation_artifacts governs "$spec_id")"
    governed_count=0
    while IFS= read -r governed_artifact; do
        [[ -n "$governed_artifact" ]] || continue
        governed_count=$((governed_count + 1))
        assert_existing_project_artifact "$governed_artifact"
        assert_contains governed_by "$governed_artifact" "$spec_id"
    done <<< "$governed_artifacts"
    if (( governed_count == 0 )); then
        printf 'Specification governs no current project artifact: %s\n' \
            "$spec_id" >&2
        exit 1
    fi

    governing_specs="$(relation_artifacts governed_by "$spec_id")"
    while IFS= read -r governing_spec; do
        [[ -n "$governing_spec" ]] || continue
        assert_existing_project_artifact "$governing_spec"
        assert_contains governs "$governing_spec" "$spec_id"
    done <<< "$governing_specs"
done

workflow_id=@/workflows/verify-specifications.donna.md
assert_exact governed_by "$workflow_id" \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
assert_contains governs @/specs/behavior/files_relations.md "$workflow_id"
assert_contains governs @/specs/general/workflows.md "$workflow_id"

spec_path_count="$(count_arguments "${spec_paths[@]}")"
printf 'Validated meta-specification governance for %d specifications, complete AGENTS.md governance, and forward/reverse governance for this workflow.\n' \
    "$spec_path_count"
```

## Repair Governance

```toml donna
id = "repair_governance"
kind = "donna.lib.request_action"
```

The specification-governance gate failed.

Standard output:

```text
{{ donna.lib.task_variable("governance_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("governance_stderr") }}
```

Query the affected specification, governed artifact, and reverse relation separately. Repair only scalable governance rules, current specification paths, the specification index, affected agent guidance, or an incorrect verifier expectation. Preserve bidirectional results and do not add a workflow-specific rule when the existing permanent-workflow pattern should apply.

After repair, `{{ donna.lib.goto("verify_index_and_structure") }}`. If ownership is ambiguous or conflicts with a specification, `{{ donna.lib.goto("blocked") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "capture_review_snapshot"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/verify-specifications.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"
artifact_entry_count="$(
    rg -Fxc -- '- Artifact: `./workflows/verify-specifications.donna.md`' \
        specs/general/workflows.md
)"
if (( artifact_entry_count != 1 )); then
    printf 'Workflow catalog must contain exactly one artifact entry; found %d.\n' \
        "$artifact_entry_count" >&2
    exit 1
fi

catalog_state="$(
    awk '
        /^### Verify specifications$/ {
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
        printf 'Unexpected workflow catalog state: %s\n' "$catalog_state" >&2
        exit 1
        ;;
esac

if ! rg -Fq -- '@/workflows/verify-specifications.donna.md' AGENTS.md; then
    printf 'AGENTS.md must identify the implemented specification verifier.\n' >&2
    exit 1
fi

if governance_output="$(depmesh -p llm dependencies "$workflow_id")"; then
    :
else
    status=$?
    printf 'Workflow governance query failed with exit code %d.\n' "$status" >&2
    exit "$status"
fi
for governing_spec in \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
do
    if rg -Fqx -- "- ${governing_spec}" <<< "$governance_output" >/dev/null; then
        :
    else
        status=$?
        printf 'Workflow is missing governing specification: %s\n' \
            "$governing_spec" >&2
        if (( status == 1 )); then
            exit 1
        fi
        exit "$status"
    fi
done

safe_project_diff_paths=(
    .
    ':(top,exclude)man/**'
)
git --no-optional-locks diff --check -- "${safe_project_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_project_diff_paths[@]}"

text_paths=("$workflow_path")
spec_paths_output="$(
    find specs -type f -name '*.md' -print |
        LC_ALL=C sort
)"
if [[ -n "$spec_paths_output" ]]; then
    while IFS= read -r spec_path; do
        text_paths+=("$spec_path")
    done <<< "$spec_paths_output"
else
    printf 'No Markdown specifications were found for workflow-artifact verification.\n' >&2
    exit 1
fi
for text_path in "${text_paths[@]}"; do
    if trailing_whitespace="$(
        rg -n '[[:blank:]]+$' "$text_path"
    )"; then
        printf '%s contains trailing whitespace:\n%s\n' \
            "$text_path" "$trailing_whitespace" >&2
        exit 1
    else
        status=$?
        if (( status != 1 )); then
            printf 'Trailing-whitespace scan failed for %s with exit code %d.\n' \
                "$text_path" "$status" >&2
            exit "$status"
        fi
    fi
    if crlf_violations="$(
        LC_ALL=C rg -n $'\r$' "$text_path"
    )"; then
        printf '%s contains CRLF line endings:\n%s\n' \
            "$text_path" "$crlf_violations" >&2
        exit 1
    else
        status=$?
        if (( status != 1 )); then
            printf 'CRLF scan failed for %s with exit code %d.\n' \
                "$text_path" "$status" >&2
            exit "$status"
        fi
    fi
    if final_byte="$(tail -c 1 -- "$text_path")"; then
        :
    else
        status=$?
        printf 'Final-newline scan failed for %s with exit code %d.\n' \
            "$text_path" "$status" >&2
        exit "$status"
    fi
    if [[ -n "$final_byte" ]]; then
        printf '%s must end with a newline.\n' "$text_path" >&2
        exit 1
    fi
done

printf 'Workflow source, catalog lifecycle state, agent guidance, governance, and project-wide non-man diff hygiene are synchronized; Donna validation and view rendering remain agent-side checks outside an executing Donna process.\n'
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

Repair this workflow, its catalog entry, affected agent guidance, its scalable governance relation, or an in-scope formatting problem. If the project-wide non-man diff check identifies an unrelated user change, preserve it and report the blocker instead of repairing it. Keep the catalog state at `implementation in progress` until the first representative execution finishes, and validate every control-flow change before continuing.

After repair, `{{ donna.lib.goto("verify_index_and_structure") }}`. If repair requires a new contract or developer direction, `{{ donna.lib.goto("blocked") }}`.

## Capture Review Snapshot

```toml donna
id = "capture_review_snapshot"
kind = "donna.lib.run_script"
save_stdout_to = "review_snapshot"
save_stderr_to = "review_snapshot_stderr"
goto_on_success = "review_complete_result"
goto_on_failure = "repair_review_snapshot"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

emit_worktree_manifest() {
    git --no-optional-locks ls-files -z --cached --others --exclude-standard -- "$@" |
        LC_ALL=C sort -z |
        while IFS= read -r -d '' path; do
            printf '%s\0' "$path" || exit $?
            if [[ -L "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect symlink mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'symlink\0%s\0' "$mode" || exit $?
                if readlink -z -- "$path" |
                    sha256sum |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint symlink target: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -f "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect file mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'file\0%s\0' "$mode" || exit $?
                if sha256sum < "$path" |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint file contents: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -e "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect artifact mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'other\0%s\0-\0' "$mode" || exit $?
            else
                printf 'missing\0-\0-\0' || exit $?
            fi
        done
}

review_pathspecs=(
    AGENTS.md
    depmesh.toml
    specs
    workflows/verify-specifications.donna.md
)

review_fingerprint="$(
    {
        git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
            "${review_pathspecs[@]}" &&
            printf '\0review-worktree-manifest-v1\0' &&
            emit_worktree_manifest "${review_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"
if [[ ! "$review_fingerprint" =~ ^[0-9a-f]{64}$ ]]; then
    printf 'Review snapshot did not produce the expected SHA-256 identifier.\n' >&2
    exit 1
fi

printf 'v1:%s\n' "$review_fingerprint"
```

## Repair Review Snapshot

```toml donna
id = "repair_review_snapshot"
kind = "donna.lib.request_action"
```

The workflow could not capture the exact authorized edit-set state for final review.

Standard error:

```text
{{ donna.lib.task_variable("review_snapshot_stderr") }}
```

1. Diagnose only the snapshot command, required local command availability, repository state, or an incorrect verifier expectation.
2. Do not change HEAD, the Git index, package implementation, unrelated files, or any path under `man/`.
3. If a verifier repair is sufficient, `{{ donna.lib.goto("verify_index_and_structure") }}`.
4. If the review state cannot be captured without developer direction, `{{ donna.lib.goto("blocked") }}`.

## Verify Scope Stability

```toml donna
id = "verify_scope_stability"
kind = "donna.lib.run_script"
save_stdout_to = "scope_stability_stdout"
save_stderr_to = "scope_stability_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_scope_stability"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

emit_worktree_manifest() {
    git --no-optional-locks ls-files -z --cached --others --exclude-standard -- "$@" |
        LC_ALL=C sort -z |
        while IFS= read -r -d '' path; do
            printf '%s\0' "$path" || exit $?
            if [[ -L "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect symlink mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'symlink\0%s\0' "$mode" || exit $?
                if readlink -z -- "$path" |
                    sha256sum |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint symlink target: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -f "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect file mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'file\0%s\0' "$mode" || exit $?
                if sha256sum < "$path" |
                    awk '{ printf "%s%c", $1, 0 }'; then
                    :
                else
                    status=$?
                    printf 'Could not fingerprint file contents: %s\n' "$path" >&2
                    exit "$status"
                fi
            elif [[ -e "$path" ]]; then
                if mode="$(stat -c '%f' -- "$path")"; then
                    :
                else
                    status=$?
                    printf 'Could not inspect artifact mode: %s\n' "$path" >&2
                    exit "$status"
                fi
                printf 'other\0%s\0-\0' "$mode" || exit $?
            else
                printf 'missing\0-\0-\0' || exit $?
            fi
        done
}

emit_index_manifest() {
    printf 'index-stage-v3\0' || return
    git --no-optional-locks ls-files --stage -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-assume-skip-v1\0' || return
    git --no-optional-locks ls-files --stage -v -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-fsmonitor-v1\0' || return
    git --no-optional-locks ls-files --stage -f -z -- "$@" |
        LC_ALL=C sort -z || return
    printf '\0index-entry-flags-v1\0' || return
    LC_ALL=C git --no-optional-locks ls-files --stage --debug -z -- "$@" || return
    printf '\0index-resolve-undo-v1\0' || return
    git --no-optional-locks ls-files --resolve-undo -z -- "$@" |
        LC_ALL=C sort -z
}

baseline="$(
    cat <<'DONNA_SCOPE_BASELINE'
{{ donna.lib.task_variable("scope_baseline") }}
DONNA_SCOPE_BASELINE
)"
if [[ ! "$baseline" =~ ^v3:[0-9a-f]{40,64}:[0-9a-f]{64}:[0-9a-f]{64}:[0-9a-f]{64}$ ]]; then
    printf 'Saved scope baseline is missing or malformed.\n' >&2
    exit 1
fi

review_snapshot="$(
    cat <<'DONNA_REVIEW_SNAPSHOT'
{{ donna.lib.task_variable("review_snapshot") }}
DONNA_REVIEW_SNAPSHOT
)"
if [[ ! "$review_snapshot" =~ ^v1:[0-9a-f]{64}$ ]]; then
    printf 'Saved review snapshot is missing or malformed.\n' >&2
    exit 1
fi

IFS=':' read -r baseline_version baseline_head baseline_index baseline_protected baseline_initial_review \
    <<< "$baseline"
if [[ "$baseline_version" != "v3" ]]; then
    printf 'Unsupported scope-baseline version: %s\n' "$baseline_version" >&2
    exit 1
fi
IFS=':' read -r review_version baseline_review <<< "$review_snapshot"
if [[ "$review_version" != "v1" ]]; then
    printf 'Unsupported review-snapshot version: %s\n' "$review_version" >&2
    exit 1
fi

protected_pathspecs=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)AGENTS.md'
    ':(top,exclude)depmesh.toml'
    ':(top,exclude)specs/**'
    ':(top,exclude)workflows/verify-specifications.donna.md'
)

current_index="$(
    emit_index_manifest \
        . \
        ':(top,exclude)man/**' |
        sha256sum |
        awk '{ print $1 }'
)"
index_scope_description="non-man staged entries, assume-unchanged, skip-worktree, fsmonitor-valid, intent-to-add, and resolve-undo state"

current_protected="$(
    {
        git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
            "${protected_pathspecs[@]}" &&
            printf '\0protected-worktree-manifest-v1\0' &&
            emit_worktree_manifest "${protected_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"

review_pathspecs=(
    AGENTS.md
    depmesh.toml
    specs
    workflows/verify-specifications.donna.md
)
current_review="$(
    {
        git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
            "${review_pathspecs[@]}" &&
            printf '\0review-worktree-manifest-v1\0' &&
            emit_worktree_manifest "${review_pathspecs[@]}"
    } |
        sha256sum |
        awk '{ print $1 }'
)"

current_head="$(git --no-optional-locks rev-parse --verify 'HEAD^{commit}')"
mismatch=0
if [[ "$current_head" != "$baseline_head" ]]; then
    printf 'HEAD changed after baseline capture.\nBaseline: %s\nCurrent:  %s\n' \
        "$baseline_head" "$current_head" >&2
    mismatch=1
fi
if [[ "$current_index" != "$baseline_index" ]]; then
    printf 'The %s changed after baseline capture.\nBaseline: %s\nCurrent:  %s\n' \
        "$index_scope_description" "$baseline_index" "$current_index" >&2
    mismatch=1
fi
if [[ "$current_protected" != "$baseline_protected" ]]; then
    printf 'A protected tracked or non-ignored untracked non-man worktree artifact changed after baseline capture.\nBaseline: %s\nCurrent:  %s\n' \
        "$baseline_protected" "$current_protected" >&2
    mismatch=1
fi
if [[ "$current_review" != "$baseline_review" ]]; then
    printf 'The authorized specification-verification edit set changed during final review.\nBaseline: %s\nCurrent:  %s\n' \
        "$baseline_review" "$current_review" >&2
    mismatch=1
fi
if (( mismatch )); then
    exit 1
fi

if [[ "$baseline_review" == "$baseline_initial_review" ]]; then
    review_change_summary="the authorized edit set remained identical to its initial state"
else
    review_change_summary="the authorized edit set differs from its initial state and exactly matches the reviewed state"
fi
printf 'HEAD, %s, and protected tracked and non-ignored untracked non-man artifacts match the initial baseline; %s.\n' \
    "$index_scope_description" "$review_change_summary"
```

## Repair Scope Stability

```toml donna
id = "repair_scope_stability"
kind = "donna.lib.request_action"
```

The protected-scope comparison failed.

Initial baseline:

```text
{{ donna.lib.task_variable("scope_baseline") }}
```

Final-review snapshot:

```text
{{ donna.lib.task_variable("review_snapshot") }}
```

Standard output:

```text
{{ donna.lib.task_variable("scope_stability_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("scope_stability_stderr") }}
```

1. Inspect only non-`man/` status and diffs with `git --no-optional-locks`, without modifying HEAD or the Git index.
2. If HEAD, the index, a protected artifact, or the authorized edit set changed outside this workflow's reviewed repair, preserve the change and `{{ donna.lib.goto("blocked") }}` so the developer can decide how to proceed.
3. If the protected state again matches because a concurrent transient change was reverted, `{{ donna.lib.goto("verify_index_and_structure") }}` so preflight and every verification gate run again on the restored state.
4. If the comparison itself is defective, query this workflow's relations, repair only the verifier, validate the control-flow change, and `{{ donna.lib.goto("verify_index_and_structure") }}`.

## Review Complete Result

```toml donna
id = "review_complete_result"
kind = "donna.lib.request_action"
```

All deterministic specification gates passed.

Initial protected-scope baseline:

```text
{{ donna.lib.task_variable("scope_baseline") }}
```

Authorized edit-set snapshot for this review:

```text
{{ donna.lib.task_variable("review_snapshot") }}
```

1. Inspect the exact in-scope diff with `git --no-optional-locks` and the complete workflow source without modifying the Git index. If the workflow is untracked, read it directly because Git diff does not include it.
2. Reread every changed specification and enough unchanged context to confirm the complete specification set remains coherent.
3. Confirm the index contains every specification and directory exactly once, opening sections are meaningful, path conventions are correct, normative statements are clear and non-conflicting, and governance results match current paths.
4. Confirm the current HEAD, non-`man/` staged entries, assume-unchanged, skip-worktree, fsmonitor-valid, intent-to-add, and resolve-undo state, tracked or non-ignored untracked non-`man/` worktree artifacts outside the authorized edit set, and the authorized edit set itself still appear consistent with the displayed baselines. The transition after this review will enforce all comparisons again without another agent pause.
5. Confirm no package implementation behavior, generated manual, unrelated user change inside the authorized edit set, network resource, or Donna session reset entered the workflow scope.
6. While Donna is awaiting this action and no workflow operation is executing, run `donna -p llm validate @/workflows/verify-specifications.donna.md` and `donna -p llm render @/workflows/verify-specifications.donna.md --mode view`; inspect the rendered instructions and transitions.
7. Run `donna -p llm list` and confirm the workflow is discoverable with the intended title and summary.
8. Confirm the catalog, agent guidance, scalable Depmesh rules, and workflow artifact form one coherent implementation change.
9. If the complete specification contract is satisfied, `{{ donna.lib.goto("verify_scope_stability") }}` so the protected-scope comparison runs immediately before completion.
10. If a repairable discrepancy remains, `{{ donna.lib.goto("repair_review_findings") }}`.
11. If a contract conflict or ownership decision requires developer input, record an unresolved finding when required and `{{ donna.lib.goto("blocked") }}`.

## Repair Review Findings

```toml donna
id = "repair_review_findings"
kind = "donna.lib.request_action"
```

Query all configured relations separately for every affected artifact, then repair the specification, index, governance, workflow, catalog, guidance, or verifier discrepancy. Preserve unrelated changes, keep the Git index untouched, do not change package implementation behavior, and do not inspect generated manuals.

After repair, `{{ donna.lib.goto("verify_index_and_structure") }}`. If the finding cannot be resolved within the current specification contract, record it in an open-questions report when required and `{{ donna.lib.goto("blocked") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The complete specification set passed deterministic structure, index, path-style, governance, workflow-artifact, protected-scope stability, and project-wide non-man diff checks plus source-level semantic review. Report specifications repaired, index parity, normative review findings, governance evidence, protected-scope evidence, workflow validation and rendering, and any intentionally unresolved questions.

If this was the first representative execution and the catalog still says `implementation in progress`, mark it `implemented` only after this finish. Then rerun this workflow in the final catalog state and repeat specific workflow validation, view rendering, governance queries, and project-wide non-man diff checks before reporting completion.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The specification-verification workflow cannot continue without developer input. Report completed gates, the exact structural, semantic, governance, or ownership conflict, affected artifacts, any required open-questions report entry, and the decision needed. Do not inspect generated manuals, change package behavior, stage files, create commits, use the network, or reset the Donna session.
