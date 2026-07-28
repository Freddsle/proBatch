# Migrate BEC E. coli Core History

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "verify_inputs"
```

Analyze the BEC E. coli branch and the stopped split implementation, preserve one ordered manifest position for every selected source commit, run numbered workflows only for core-relevant positions with two maintainer review pauses each, record companion-only positions as deterministic reference skips, and finish with a structured residual split review.

## Verify Inputs

```toml donna
id = "verify_inputs"
kind = "donna.lib.run_script"
save_stdout_to = "input_stdout"
save_stderr_to = "input_stderr"
goto_on_success = "confirm_history_scope"
goto_on_failure = "repair_inputs"
timeout = 120
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

source_repo="/home/yuliya/repos/other/proBatch"
split_repo="/home/yuliya/repos/other/proBatch-core-split"
api_map="/home/yuliya/repos/other/proBatch/archive/stopped-probatch-split-2026-07-27/API_FUNCTIONS.md"
base="ba6ee246eace090e71baa7aba302ca64e76ddb32"
sync_merge="5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab"
tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"
split_implementation="29a7478dc7deea846a2c1ff1abd25a881e6f87db"
split_stop_note="a04c9a29a3a4ba9f719d3b9c778616f3dd77903b"
split_preflight="49cee7cc978fbb149c262a5a783face32dd1d135"
expected_api_hash="21649f0a9b7eb178886ed0e466f4dcc1a757b0ba7ac5c32f9863cd82c82c2278"
expected_ancestry_hash="788f5593e198e01a3292b5f82bd76bc123e3a48f191fa173b62571f7ee62574c"

test -f "$api_map"
test "$(git -C "$source_repo" rev-parse --is-inside-work-tree)" = "true"
test "$(git -C "$split_repo" rev-parse --is-inside-work-tree)" = "true"

test "$(sha256sum "$api_map" | awk '{print $1}')" = "$expected_api_hash"
test "$(git -C "$source_repo" rev-parse refs/heads/BEC_ecoli_data)" = "$tip"
git -C "$source_repo" cat-file -e "$base^{commit}"
git -C "$source_repo" cat-file -e "$sync_merge^{commit}"
git -C "$source_repo" cat-file -e "$tip^{commit}"
git -C "$split_repo" cat-file -e "$split_implementation^{commit}"
git -C "$split_repo" cat-file -e "$split_stop_note^{commit}"
git -C "$split_repo" cat-file -e "$split_preflight^{commit}"
git -C "$source_repo" merge-base --is-ancestor "$base" "$tip"
test "$(git -C "$source_repo" rev-parse "$sync_merge^2")" = "$base"
test "$(git -C "$split_repo" rev-parse refs/heads/split/probatch-core)" = "$split_stop_note"
test "$(git -C "$split_repo" rev-parse "$split_stop_note^")" = "$split_implementation"
test "$(git -C "$split_repo" rev-parse "$split_implementation^")" = "$split_preflight"
test "$(git -C "$split_repo" rev-parse "$split_preflight^")" = "$base"
test "$(git -C "$source_repo" rev-list --reverse --ancestry-path "$base..$tip" | sha256sum | awk '{print $1}')" = "$expected_ancestry_hash"

full_difference_count="$(git -C "$source_repo" rev-list --count "$base..$tip")"
descendant_path_count="$(git -C "$source_repo" rev-list --count --ancestry-path "$base..$tip")"
post_sync_count="$(git -C "$source_repo" rev-list --count "$sync_merge..$tip")"
pre_sync_merge_base="$(git -C "$source_repo" merge-base "$sync_merge^1" "$base")"
pre_sync_branch_count="$(git -C "$source_repo" rev-list --count "$pre_sync_merge_base..$sync_merge^1")"

test "$full_difference_count" -eq 312
test "$descendant_path_count" -eq 30
test "$post_sync_count" -eq 29
test "$pre_sync_branch_count" -eq 282

printf 'source branch: refs/heads/BEC_ecoli_data at %s\n' "$tip"
printf 'source baseline: %s\n' "$base"
printf 'synchronizing merge: %s\n' "$sync_merge"
printf 'stopped split implementation: %s\n' "$split_implementation"
printf 'stopped split provenance note: %s\n' "$split_stop_note"
printf 'stopped split preflight: %s\n' "$split_preflight"
printf 'all commits reachable from tip but not base: %s\n' "$full_difference_count"
printf 'commits on a descendant path from base: %s\n' "$descendant_path_count"
printf 'linear commits after the synchronizing merge: %s\n' "$post_sync_count"
printf 'older divergent branch commits introduced before synchronization: %s\n' "$pre_sync_branch_count"
printf 'source worktree status is intentionally ignored; historical source reads must use pinned Git objects\n'
```

## Repair Inputs

```toml donna
id = "repair_inputs"
kind = "donna.lib.request_action"
```

The pinned source verification failed.

Standard output:

```text
{{ donna.lib.task_variable("input_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("input_stderr") }}
```

1. Diagnose changed or missing repositories, refs, commits, or the API map without modifying either reference repository.
2. Do not substitute a new branch tip, base, synchronizing merge, API map, split preflight, split implementation, split stop-note, recorded hash, or expected history count without the developer's explicit approval.
3. If the pinned inputs have been restored or explicitly revised, `{{ donna.lib.goto("verify_inputs") }}`.
4. If source identity cannot be resolved safely, `{{ donna.lib.goto("blocked") }}`.

## Confirm History Scope

```toml donna
id = "confirm_history_scope"
kind = "donna.lib.request_action"
```

The source topology is:

```text
{{ donna.lib.task_variable("input_stdout") }}
```

The requested name `BEC_ecoli_databranch` does not exist. The exact local branch is `refs/heads/BEC_ecoli_data` at pinned tip `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.

The base `ba6ee246eace090e71baa7aba302ca64e76ddb32` entered the BEC branch as the second parent of merge `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`. Consequently, the raw Git set `base..tip` contains 312 commits: 282 older divergent branch commits, the synchronizing merge, and 29 commits made after synchronization. It is not a linear list of 312 commits authored after the named base.

Choose the migration universe explicitly:

1. Recommended: create one foundation workflow for the core-owned API present at the synchronizing merge, then preserve one ordered position for each of the 29 linear commits after that merge. Materialize executable children only for core-relevant positions, while retaining unambiguously companion-only commits as reference skips. This preserves per-commit review where core work exists while still accounting for the complete source history. With this scope, `{{ donna.lib.goto("analyze_post_sync_history") }}`.
2. Literal Git difference: create one workflow for every one of the 312 commits reachable from the tip but not the base, including the 282 older divergent commits and merge histories. With this scope, `{{ donna.lib.goto("analyze_full_difference") }}`.
3. If neither scope represents the intended history, `{{ donna.lib.goto("blocked") }}` and report the exact alternative ref/range or ordering rule required.

Do not complete this action request until the developer explicitly chooses a scope.

## Analyze Post-Synchronization History

```toml donna
id = "analyze_post_sync_history"
kind = "donna.lib.request_action"
```

1. Read the complete API ownership map at `/home/yuliya/repos/other/proBatch/archive/stopped-probatch-split-2026-07-27/API_FUNCTIONS.md`.
2. Treat both reference repositories as immutable. Inspect pinned trees and diffs only with `git -C <repo> show`, `diff`, `ls-tree`, or `grep` against full commit IDs, or with a read-only archive created under `/tmp`. Never inspect an external working-tree version and never run checkout, switch, reset, stash, clean, or any command that writes to either reference repository.
3. Inspect the non-generated implementation at pinned split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` and read `SPLIT_ATTEMPT.md` from pinned stop-note commit `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b` through Git object commands in `/home/yuliya/repos/other/proBatch-core-split`. Treat the note and implementation as rejected historical evidence, not as an accepted patch. Inventory the complete non-generated split change from baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32` through the mode-only preflight `49cee7cc978fbb149c262a5a783face32dd1d135` to the implementation commit, including DESCRIPTION, R and Roxygen sources, tests, vignettes, data, `inst/` resources, and project configuration.
4. Do not read, enumerate, compare, or transfer generated `man/` files in either reference repository, and do not transfer or validate `NAMESPACE`.
5. Define a foundation snapshot at merge `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab` that covers every API function assigned to core, its editable Roxygen source, focused tests, and required package metadata without importing companion-owned providers or benchmark code.
6. Analyze the exact output of `git rev-list --reverse --topo-order 5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab..e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`, parent before child and oldest to newest. It must contain 29 commits, begin with `42544a21f10ca6960d3e4c44d2833f764054d721`, and end with `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`. For every commit record the full SHA, parent, date, complete subject and body, non-generated changed files, relevant hunks and intent, core/companion/mixed classification, core functions affected, focused tests, package dependencies, and stopped-split comparison targets.
7. At hunk level, separate core behavior from all 26 companion-owned exports. This includes provider-specific PRONE, BERT, NormAE, PLSDA, RUVIIIC, mComBat, omicsGMF, and missForest adapters plus `calculate_classification_metrics`, `calculate_variance_partition`, `prepare_variance_partition_df`, `plot_variance_partition`, `plot_variance_partition.df`, `plot_intragroup_variation`, `plot_TSNE`, and `plot_UMAP`. Retain independent core hunks from mixed files such as `R/proteome_wide_diagnostics.R`.
8. Account explicitly for the load-order ambiguity from `correct_batch_effects_old.R`; determine effective behavior rather than choosing by filename recency.
9. Write the persistent ordered manifest to `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`. Include one foundation entry, 29 commit entries, and a final residual-review entry. Each entry must state why it is executable, reference-only, skipped, or split. Add the exact machine-readable header `<!-- donna-migration scope=post-sync last-prefix=032 entries=31 workflows=24 commits=29 references=7 -->`. Add one machine-readable entry per allocated child position, in numeric order, using `<!-- donna-entry kind=<foundation|commit|residual> id=<unique-id> slot=<NNN> mode=<workflow|reference-only> workflow=<artifact-id-or-dash> sha=<full-commit-or-dash> -->`. Use `kind=commit` and the full source SHA for exactly the 29 commit entries. Slots `005`, `018`, `024`, `025`, `026`, `029`, and `030` must use `mode=reference-only workflow=-`; all other positions must use `mode=workflow`, including the generated-only core-documentation entry at slot `017`.
10. When the audit is complete and every source commit is accounted for exactly once, `{{ donna.lib.goto("materialize_workflows") }}`.
11. If any commit, ownership decision, or source identity remains ambiguous, `{{ donna.lib.goto("blocked") }}`.

## Analyze Full Git Difference

```toml donna
id = "analyze_full_difference"
kind = "donna.lib.request_action"
```

1. Read the complete API ownership map and inspect the pinned BEC implementation, stopped-split implementation, and `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` under the same immutable Git-object, complete non-generated artifact-family, and generated-file restrictions as the post-synchronization path.
2. Analyze all 312 commits in `ba6ee246eace090e71baa7aba302ca64e76ddb32..e2bb18547c73f1c471fc1afcb3facbd8bea5fa92` using `--reverse --topo-order`. Do not sort by author date.
3. Preserve the DAG: record every parent and, for each merge, the exact parent comparison used to identify the transferable effect. Never treat a merge as an ordinary single-parent patch.
4. For every commit record the full SHA, parents, date, subject and body, non-generated changed files, relevant hunks and intent, core/companion/mixed classification, core functions affected, tests, dependencies, and stopped-split comparison targets.
5. Exclude all 26 companion-owned exports named in the post-synchronization analysis at hunk level while retaining independent core hunks from mixed commits. Do not cherry-pick commits wholesale.
6. Account explicitly for the effective installed behavior and load-order ambiguity from `correct_batch_effects_old.R`.
7. Write the persistent ordered manifest to `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}` with all 312 commits and one final residual-review entry. Explain the chosen topological order and merge comparators. Add the exact machine-readable header `<!-- donna-migration scope=full-difference last-prefix=314 entries=313 workflows=313 commits=312 references=0 -->` and one machine-readable entry in the same slot/mode format required by the post-synchronization path. Use `mode=workflow`, `kind=commit`, and the full source SHA for exactly the 312 commit entries.
8. When every commit is accounted for exactly once, `{{ donna.lib.goto("materialize_workflows") }}`.
9. If the DAG cannot be linearized without changing the intended review semantics, `{{ donna.lib.goto("blocked") }}`.

## Materialize Numbered Workflows

```toml donna
id = "materialize_workflows"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}` and query Depmesh for every artifact that will be edited.
2. Allocate zero-padded three-digit slots in manifest order after this `001_` parent. Use a unique 12-hexadecimal-character commit abbreviation and filenames of the form `NNN_<12-char-sha>_<concise-slug>.donna.md` for executable commit workflows. If the chosen manifest has a foundation entry, use the next slot for `NNN_foundation_core_api.donna.md`. Give the final residual workflow the last allocated slot. Never renumber or reuse an allocated slot, including one reserved by a reference-only entry.
3. Create permanent child workflows only for entries with `mode=workflow`. For the selected post-synchronization manifest, slots `005`, `018`, `024`, `025`, `026`, `029`, and `030` are exact companion-only `mode=reference-only workflow=-` records and must have no artifacts, catalog entries, or maintainer pauses. Keep slot `017` executable as the explicit generated-only core-documentation exception. A reference-only commit retains its full pinned evidence and concise disposition in the manifest and receives only a deterministic `reference-only` progress outcome when reached.
4. Every foundation or commit child must implement this control flow with stable operation IDs:
   - inspect the exact pinned commit/tree and recheck core ownership;
   - query all Depmesh relations before editing destination artifacts;
   - port only the BEC core-owned source behavior and focused tests, or record why no destination change is warranted;
   - run deterministic source-port verification with a repair loop;
   - pause at `review_bec_port` for manual developer analysis and commit, and forbid completing the action request until the developer explicitly resumes;
   - after resume, run `reverify_bec_port` so any maintainer edits are checked again before the split comparison, routing failures back through repair;
   - compare the same behavior with pinned split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db`;
   - apply only justified split hardening, or record why no split-stage change is warranted;
   - run deterministic split-stage verification with a repair loop;
   - pause at `review_split_adjustment` for manual developer analysis and commit, again requiring explicit resume;
   - after resume, run `reverify_split_adjustment` so maintainer edits are checked before finish, routing failures back through repair;
   - finish with an outcome that the parent can record.
5. Child workflows must be safe to resume and rerun: detect already-equivalent code, preserve unrelated user changes, never cherry-pick, never create or amend a commit, never draft a commit message unless asked, and never read, edit, generate, lint, validate, enumerate, or diff-review `man/`. Make documentation changes only in Roxygen2 comments under `R/`; leave `NAMESPACE` and documentation generation to the maintainer. Treat both external repositories as immutable and require every source and split comparison to use pinned Git object commands or a read-only `/tmp` archive. Forbid checkout, switch, reset, stash, clean, and all writes in the external repositories.
6. Each executable child must name its exact source SHA and comparator, relevant core functions and hunks, expected destination files, focused tests, stopped-split files, dependency questions, and known skip rationale from the manifest.
7. Create the final residual workflow with a deterministic inventory check and repair loop followed by an explicit `review_residual_reports` action that remains pending until the developer resumes. It must read the pinned `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md`, account for the overall non-generated `ba6ee246eace090e71baa7aba302ca64e76ddb32..29a7478dc7deea846a2c1ff1abd25a881e6f87db` split including preflight `49cee7cc978fbb149c262a5a783face32dd1d135`, and use pinned Git objects to build an exhaustive editable/reference path inventory for `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92..29a7478dc7deea846a2c1ff1abd25a881e6f87db`, excluding `man/**` and generated `NAMESPACE`. Its deterministic inventory operation must use ID `verify_residual_inventory`. The immutable delta contains 99 paths: 11 added, 38 deleted, and 50 modified, with name-status SHA-256 `28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38`. Account for every path exactly once, compare it with the completed core, and classify it as `required`, `recommended`, `equivalent`, `excluded`, or `decision`. In the residual review report, include exactly one marker per path, in Git name-status order, using `<!-- split-path status=<A|D|M> path=<path> class=<classification> -->`; the extracted status/path lines must reproduce the pinned count and checksum. Add exactly one `<!-- bench-api symbol=<exact-name> disposition=companion -->` marker for each of the 26 companion-owned exports. Explicitly assess the split-only `R/identifiers.R`, `R/matrix_adapter.R`, `R/registry.R`, and `R/step_result.R` layers; `tests/testthat/helper-source-root.R`, `test-identifiers.R`, `test-lineage.R`, `test-matrix_adapter.R`, `test-registry.R`, `test-step_result.R`, and `test-symbol-ownership.R`; the new `pb_apply_matrix_method`, `pb_step_result`, and `pb_unregister_steps` APIs; all 26 companion-owned exports; provider ownership removals; duplicate-definition consolidation; and retained/deleted test coverage. In the remaining-change plan, add exactly one `<!-- plan-path path=<path> class=<required|recommended|decision> destination=<nonspace-target> order=<positive-integer> reason=<nonspace-reason-key> -->` marker for every residual path with an actionable class, and explain each marker with a concise destination, reason, dependency/test impact, and execution order. If no actionable paths exist, include exactly `<!-- plan-empty -->`. It must create:
   - `{{ donna.lib.path("@/workflows/reports/bec-ecoli-core-residual-review.md") }}`
   - `{{ donna.lib.path("@/workflows/reports/bec-ecoli-core-remaining-change-plan.md") }}`
   The reports must be concise, logically structured, and include source symbol/path, destination, originating evidence, reason, dependencies, tests, ordering, and explicit exclusions.
8. Add every exact permanent child artifact to `{{ donna.lib.path("@/specs/general/workflows.md") }}` with state `implementation in progress`; do not catalog reference-only records as workflows. Update `{{ donna.lib.path("@/AGENTS.md") }}` for executable and reserved slots plus maintainer-controlled pauses. Existing Depmesh wildcard governance should be verified; change it only if the remaining artifacts are not discovered.
9. When all workflows and catalog entries exist, `{{ donna.lib.goto("validate_workflows") }}`.
10. If a workflow cannot be made specific enough for safe independent execution, `{{ donna.lib.goto("blocked") }}`.

## Validate Workflows

```toml donna
id = "validate_workflows"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_validation_stdout"
save_stderr_to = "workflow_validation_stderr"
goto_on_success = "review_generated_workflows"
goto_on_failure = "repair_workflows"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

manifest="workflows/bec-ecoli-core-migration-manifest.md"
catalog="specs/general/workflows.md"
source_repo="/home/yuliya/repos/other/proBatch"
base="ba6ee246eace090e71baa7aba302ca64e76ddb32"
sync_merge="5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab"
tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"

test -f "$manifest"

header_pattern='^<!-- donna-migration scope=(post-sync|full-difference) last-prefix=[0-9]{3} entries=[0-9]+ workflows=[0-9]+ commits=[0-9]+ references=[0-9]+ -->$'
mapfile -t headers < <(grep -E "$header_pattern" "$manifest")
test "$(grep -Ec "$header_pattern" "$manifest")" -eq 1
read -r scope last_prefix entry_count workflow_count commit_count reference_count < <(
    printf '%s\n' "${headers[0]}" |
        sed -E 's/^<!-- donna-migration scope=([^ ]+) last-prefix=([0-9]{3}) entries=([0-9]+) workflows=([0-9]+) commits=([0-9]+) references=([0-9]+) -->$/\1 \2 \3 \4 \5 \6/'
)

case "$scope" in
    post-sync)
        test "$last_prefix" = "032"
        test "$entry_count" -eq 31
        test "$workflow_count" -eq 24
        test "$commit_count" -eq 29
        test "$reference_count" -eq 7
        source_range="$sync_merge..$tip"
        expected_foundations=1
        ;;
    full-difference)
        test "$last_prefix" = "314"
        test "$entry_count" -eq 313
        test "$workflow_count" -eq 313
        test "$commit_count" -eq 312
        test "$reference_count" -eq 0
        source_range="$base..$tip"
        expected_foundations=0
        ;;
    *)
        exit 1
        ;;
esac

entry_pattern='^<!-- donna-entry kind=(foundation|commit|residual) id=[a-z0-9-]+ slot=[0-9]{3} mode=(workflow|reference-only) workflow=(-|@/workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md) sha=(-|[0-9a-f]{40}) -->$'
mapfile -t entries < <(grep -E "$entry_pattern" "$manifest")
test "$(grep -Ec "$entry_pattern" "$manifest")" -eq "$entry_count"
test "$(grep -c '^<!-- donna-entry ' "$manifest")" -eq "$entry_count"

entry_ids=()
entry_slots=()
entry_kinds=()
entry_modes=()
entry_workflows=()
workflow_artifacts=()
workflow_kinds=()
manifest_shas=()
reference_ids=()
reference_slots=()
reference_shas=()
foundation_count=0
residual_count=0
actual_workflow_count=0
actual_reference_count=0

for index in "${!entries[@]}"; do
    entry="${entries[$index]}"
    kind="$(printf '%s\n' "$entry" | sed -E 's/^<!-- donna-entry kind=([^ ]+) .*/\1/')"
    entry_id="$(printf '%s\n' "$entry" | sed -E 's/.* id=([^ ]+) slot=.*/\1/')"
    slot="$(printf '%s\n' "$entry" | sed -E 's/.* slot=([0-9]{3}) mode=.*/\1/')"
    mode="$(printf '%s\n' "$entry" | sed -E 's/.* mode=([^ ]+) workflow=.*/\1/')"
    artifact="$(printf '%s\n' "$entry" | sed -E 's/.* workflow=([^ ]+) sha=.*/\1/')"
    sha="$(printf '%s\n' "$entry" | sed -E 's/.* sha=([^ ]+) -->$/\1/')"
    expected_slot="$(printf '%03d' "$((index + 2))")"

    test "$slot" = "$expected_slot"
    entry_ids+=("$entry_id")
    entry_slots+=("$slot")
    entry_kinds+=("$kind")
    entry_modes+=("$mode")
    entry_workflows+=("$artifact")

    case "$kind" in
        foundation)
            foundation_count=$((foundation_count + 1))
            test "$sha" = "-"
            test "$mode" = "workflow"
            ;;
        commit)
            test "$sha" != "-"
            manifest_shas+=("$sha")
            ;;
        residual)
            residual_count=$((residual_count + 1))
            test "$sha" = "-"
            test "$mode" = "workflow"
            ;;
    esac

    case "$mode" in
        workflow)
            actual_workflow_count=$((actual_workflow_count + 1))
            test "$artifact" != "-"
            workflow="${artifact#@/}"
            test -f "$workflow"
            filename="${workflow#workflows/}"
            test "${filename%%_*}" = "$slot"
            grep -Fq "29a7478dc7deea846a2c1ff1abd25a881e6f87db" "$workflow"
            workflow_artifacts+=("$artifact")
            workflow_kinds+=("$kind")

            if test "$kind" = "commit"; then
                grep -Fq "$sha" "$workflow"
                abbreviation="${sha:0:12}"
                printf '%s\n' "$workflow" |
                    grep -Eq "^workflows/${slot}_${abbreviation}_[a-z0-9_-]+\\.donna\\.md$"
            fi
            ;;
        reference-only)
            actual_reference_count=$((actual_reference_count + 1))
            test "$kind" = "commit"
            test "$artifact" = "-"
            test "$(find workflows -maxdepth 1 -type f -name "${slot}_*.donna.md" -print | wc -l)" -eq 0
            reference_ids+=("$entry_id")
            reference_slots+=("$slot")
            reference_shas+=("$sha")
            ;;
    esac
done

test "$actual_workflow_count" -eq "$workflow_count"
test "$actual_reference_count" -eq "$reference_count"
test "$foundation_count" -eq "$expected_foundations"
test "$residual_count" -eq 1
test "$(printf '%s\n' "${entry_ids[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${entry_slots[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${workflow_artifacts[@]}" | sort -u | wc -l)" -eq "$workflow_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | sort -u | wc -l)" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | cut -c1-12 | sort -u | wc -l)" -eq "$commit_count"

if test "$scope" = "post-sync"; then
    test "${entry_kinds[0]}" = "foundation"
    test "${entry_modes[0]}" = "workflow"
    test "${entry_workflows[0]}" = "@/workflows/002_foundation_core_api.donna.md"
    grep -Fq "$sync_merge" "${entry_workflows[0]#@/}"
    grep -Fq "$base" "${entry_workflows[0]#@/}"
    for ((index = 1; index <= commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
    done
    test "${entry_kinds[$((entry_count - 1))]}" = "residual"
    test "${reference_slots[*]}" = "005 018 024 025 026 029 030"
    test "${reference_shas[*]}" = "65f70a46c4cf44e2717744aeafbf8acbe83b0378 b8a262b4256966d60e1e8452ebde7a1bf471b4af 4a06d99949114b2804b8d34492a288872fb611ed 6601232c69b44c507f2e9f63d256836e655f7973 96d38eb7449e4d38c0d6a3fffd66e17f145f669f 0b3c3b55403f2bb7342a37236dd30f3d9b3544e9 20e76c9b9b28c0ec98faa63c3f9382c1347301b9"
    test "$(grep -Fxc '<!-- donna-entry kind=commit id=post-sync-015-4540aca9-generated-missing-docs slot=017 mode=workflow workflow=@/workflows/017_4540aca9182c_generated_missing_docs.donna.md sha=4540aca9182c6708fe9bda0b8fc33d2cf8c13e57 -->' "$manifest")" -eq 1
else
    for ((index = 0; index < commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
        test "${entry_modes[$index]}" = "workflow"
    done
    test "${entry_kinds[$commit_count]}" = "residual"
fi

mapfile -t expected_shas < <(git -C "$source_repo" rev-list --reverse --topo-order "$source_range")
test "$(printf '%s\n' "${expected_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
for index in "${!expected_shas[@]}"; do
    test "${manifest_shas[$index]}" = "${expected_shas[$index]}"
done

mapfile -t numbered_files < <(find workflows -maxdepth 1 -type f -name '[0-9][0-9][0-9]_*.donna.md' -print | sort)
numbered_count="$(printf '%s\n' "${numbered_files[@]}" | wc -l)"
test "$numbered_count" -eq "$((workflow_count + 1))"
test "${numbered_files[0]}" = "workflows/001_bec-ecoli-core-migration.donna.md"
expected_numbered_files=("workflows/001_bec-ecoli-core-migration.donna.md")
for artifact in "${workflow_artifacts[@]}"; do
    expected_numbered_files+=("${artifact#@/}")
done
for index in "${!numbered_files[@]}"; do
    test "${numbered_files[$index]}" = "${expected_numbered_files[$index]}"
done

catalog_pattern='^- Artifact: `\./workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md`$'
test "$(grep -Ec "$catalog_pattern" "$catalog")" -eq "$((workflow_count + 1))"
governs_output="$(depmesh -p llm dependencies --relation governs @/specs/general/workflows.md)"
for workflow in "${numbered_files[@]}"; do
    artifact="@/${workflow}"
    catalog_line="- Artifact: \`./${workflow}\`"
    dependency_output="$(depmesh -p llm dependencies "$artifact")"

    test "$(grep -Fxc -- "$catalog_line" "$catalog")" -eq 1
    awk -v target="$catalog_line" '
        $0 == target {
            getline
            found = ($0 == "- State: implementation in progress")
        }
        END {
            exit(found ? 0 : 1)
        }
    ' "$catalog"
    printf '%s\n' "$dependency_output" | grep -Fxq -- "- @/specs/behavior/files_relations.md"
    printf '%s\n' "$dependency_output" | grep -Fxq -- "- @/specs/general/workflows.md"
    test "$(printf '%s\n' "$dependency_output" | grep -c '^- @/specs/')" -eq 2
    printf '%s\n' "$governs_output" | grep -Fxq -- "- $artifact"

    donna -p llm validate "$artifact"
    donna -p llm render "$artifact" --mode view >/dev/null
done

for index in "${!workflow_artifacts[@]}"; do
    workflow="${workflow_artifacts[$index]#@/}"
    kind="${workflow_kinds[$index]}"

    if test "$kind" = "foundation" || test "$kind" = "commit"; then
        test "$(grep -Fxc 'id = "review_bec_port"' "$workflow")" -eq 1
        test "$(grep -Fxc 'id = "reverify_bec_port"' "$workflow")" -eq 1
        test "$(grep -Fxc 'id = "review_split_adjustment"' "$workflow")" -eq 1
        test "$(grep -Fxc 'id = "reverify_split_adjustment"' "$workflow")" -eq 1
    else
        test "$(grep -Fxc 'id = "verify_residual_inventory"' "$workflow")" -eq 1
        test "$(grep -Fxc 'id = "review_residual_reports"' "$workflow")" -eq 1
    fi
done

donna -p llm validate --all
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair Workflows

```toml donna
id = "repair_workflows"
kind = "donna.lib.request_action"
```

The generated workflow or governance validation failed.

Standard output:

```text
{{ donna.lib.task_variable("workflow_validation_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("workflow_validation_stderr") }}
```

1. Repair only the manifest, generated workflows, workflow catalog, agent guidance, or affected governance rules.
2. Preserve the selected scope, sequence numbers, source SHAs, two maintainer pauses, and generated-file prohibitions.
3. When repaired, `{{ donna.lib.goto("validate_workflows") }}`.
4. If repair requires changing the selected history scope or API ownership, `{{ donna.lib.goto("blocked") }}`.

## Review Generated Workflows

```toml donna
id = "review_generated_workflows"
kind = "donna.lib.request_action"
```

The manifest, numbered child workflows, catalog entries, and agent guidance validate and render successfully.

1. Confirm and report the exact slot accounting: 31 ordered entries through prefix `032`, 24 executable children, seven reference-only entries with reserved gaps, 25 numbered workflow artifacts including the parent, unique entry IDs and commit SHAs, parent-before-child SHA order, exact manifest-to-executable-child-to-catalog mapping, exact Depmesh governance, and validation/render result.
2. Confirm every executable foundation/commit child has both maintainer review gates and both post-resume reverification gates, slot `017` remains executable, every reference-only slot has no workflow or catalog artifact, and the residual child has a deterministic 99-path inventory gate, repair loop, and explicit report-review gate.
3. Report the selected scope, classifications, generated files, catalog changes, and validation results to the developer.
4. Give the developer time to inspect and manually commit the workflow suite separately from package implementation changes.
5. Do not create or amend that commit, and do not complete this action request until the developer explicitly asks to resume.
6. If the developer accepts the workflow suite and resumes, `{{ donna.lib.goto("run_next_commit_workflow") }}`.
7. If the developer requests workflow corrections, `{{ donna.lib.goto("repair_workflows") }}`.
8. If the workflow suite is rejected or requires a different history scope, `{{ donna.lib.goto("blocked") }}`.

## Run Next Commit Workflow

```toml donna
id = "run_next_commit_workflow"
kind = "donna.lib.request_action"
```

1. Read the ordered manifest and `.session/donna/bec-ecoli-core-migration-progress.md`, creating the progress file with `apply_patch` if it does not yet exist. Record each processed position exactly once with `<!-- donna-complete id=<manifest-entry-id> slot=<NNN> workflow=<artifact-id-or-dash> outcome=<source-and-split|source-only|split-only|no-change|reference-only|reports> source-commit=<40-hex-or-none> split-commit=<40-hex-or-none> -->`. The commit fields are destination-repository commits confirmed by the developer at the corresponding manual review pauses, not the immutable reference commits. Use `none` when that stage produced no developer commit. Executable foundation and commit children may use only the first four outcomes, consistently with which review commits exist. A reference-only entry must use `workflow=- outcome=reference-only source-commit=none split-commit=none`; reserve `reports` for the residual child, with `source-commit=none` and its developer-confirmed report commit (or `none`) in `split-commit`.
2. Select exactly the first unfinished foundation or commit position in manifest order. If it has `mode=reference-only`, verify that it is one of the seven allowed slot/SHA records, add its deterministic reference-only completion marker with `apply_patch`, do not run Donna or modify package artifacts, and `{{ donna.lib.goto("run_next_commit_workflow") }}`.
3. If the selected position has `mode=workflow`, run its exact child with `donna -p llm run @/workflows/<exact-numbered-name>.donna.md`.
4. Complete an executable child before completing this parent action request. The child will stop twice for the developer to inspect and manually commit. Never complete either child review request until the developer explicitly asks to resume.
5. After the executable child finishes and both post-resume verification gates pass, record its exact workflow, SHA or foundation snapshot, source-stage result, split-stage result, verification, developer-confirmed commits or no-change decisions, and the required machine-readable completion marker in the session progress file.
6. If another foundation or commit position remains, `{{ donna.lib.goto("run_next_commit_workflow") }}`.
7. If every foundation and commit position has either an executable completion or a reference-only completion, `{{ donna.lib.goto("run_residual_workflow") }}`.
8. If a child reports an unresolved blocker, `{{ donna.lib.goto("blocked") }}`.

## Run Residual Workflow

```toml donna
id = "run_residual_workflow"
kind = "donna.lib.request_action"
```

1. Run the final numbered residual split workflow from the manifest as a child.
2. Complete the child, including its report review gate, before completing this parent action request.
3. Confirm both required report files exist, then record the residual entry with its exact slot and workflow, `outcome=reports`, `source-commit=none`, and the developer-confirmed report commit—or `split-commit=none` for an explicit no-new-commit decision—in the same machine-readable completion-marker schema.
4. If the residual review and reports are complete, `{{ donna.lib.goto("verify_completion") }}`.
5. If the residual review is blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify Completion

```toml donna
id = "verify_completion"
kind = "donna.lib.run_script"
save_stdout_to = "completion_stdout"
save_stderr_to = "completion_stderr"
goto_on_success = "finalize_catalog"
goto_on_failure = "repair_completion"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

manifest="workflows/bec-ecoli-core-migration-manifest.md"
progress=".session/donna/bec-ecoli-core-migration-progress.md"
residual_report="workflows/reports/bec-ecoli-core-residual-review.md"
remaining_plan="workflows/reports/bec-ecoli-core-remaining-change-plan.md"
source_repo="/home/yuliya/repos/other/proBatch"
base="ba6ee246eace090e71baa7aba302ca64e76ddb32"
sync_merge="5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab"
tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"
split_implementation="29a7478dc7deea846a2c1ff1abd25a881e6f87db"

test -f "$manifest"
test -s "$residual_report"
test -s "$remaining_plan"
test -f "$progress"

header_pattern='^<!-- donna-migration scope=(post-sync|full-difference) last-prefix=[0-9]{3} entries=[0-9]+ workflows=[0-9]+ commits=[0-9]+ references=[0-9]+ -->$'
mapfile -t headers < <(grep -E "$header_pattern" "$manifest")
test "$(grep -Ec "$header_pattern" "$manifest")" -eq 1
read -r scope last_prefix entry_count workflow_count commit_count reference_count < <(
    printf '%s\n' "${headers[0]}" |
        sed -E 's/^<!-- donna-migration scope=([^ ]+) last-prefix=([0-9]{3}) entries=([0-9]+) workflows=([0-9]+) commits=([0-9]+) references=([0-9]+) -->$/\1 \2 \3 \4 \5 \6/'
)

case "$scope" in
    post-sync)
        test "$last_prefix" = "032"
        test "$entry_count" -eq 31
        test "$workflow_count" -eq 24
        test "$commit_count" -eq 29
        test "$reference_count" -eq 7
        source_range="$sync_merge..$tip"
        expected_foundations=1
        ;;
    full-difference)
        test "$last_prefix" = "314"
        test "$entry_count" -eq 313
        test "$workflow_count" -eq 313
        test "$commit_count" -eq 312
        test "$reference_count" -eq 0
        source_range="$base..$tip"
        expected_foundations=0
        ;;
    *)
        exit 1
        ;;
esac

entry_pattern='^<!-- donna-entry kind=(foundation|commit|residual) id=[a-z0-9-]+ slot=[0-9]{3} mode=(workflow|reference-only) workflow=(-|@/workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md) sha=(-|[0-9a-f]{40}) -->$'
completion_pattern='^<!-- donna-complete id=[a-z0-9-]+ slot=[0-9]{3} workflow=(-|@/workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md) outcome=(source-and-split|source-only|split-only|no-change|reference-only|reports) source-commit=(none|[0-9a-f]{40}) split-commit=(none|[0-9a-f]{40}) -->$'
mapfile -t entries < <(grep -E "$entry_pattern" "$manifest")
mapfile -t completions < <(grep -E "$completion_pattern" "$progress")
test "$(grep -Ec "$entry_pattern" "$manifest")" -eq "$entry_count"
test "$(grep -c '^<!-- donna-entry ' "$manifest")" -eq "$entry_count"
test "$(grep -Ec "$completion_pattern" "$progress")" -eq "$entry_count"
test "$(grep -c '^<!-- donna-complete ' "$progress")" -eq "$entry_count"

entry_ids=()
entry_slots=()
entry_kinds=()
entry_modes=()
entry_workflows=()
workflow_artifacts=()
manifest_shas=()
reference_slots=()
reference_shas=()
foundation_count=0
residual_count=0
actual_workflow_count=0
actual_reference_count=0
review_commits=()
review_commit_count=0

for index in "${!entries[@]}"; do
    entry="${entries[$index]}"
    completion="${completions[$index]}"
    entry_kind="$(printf '%s\n' "$entry" | sed -E 's/^<!-- donna-entry kind=([^ ]+) .*/\1/')"
    entry_id="$(printf '%s\n' "$entry" | sed -E 's/.* id=([^ ]+) slot=.*/\1/')"
    entry_slot="$(printf '%s\n' "$entry" | sed -E 's/.* slot=([0-9]{3}) mode=.*/\1/')"
    entry_mode="$(printf '%s\n' "$entry" | sed -E 's/.* mode=([^ ]+) workflow=.*/\1/')"
    entry_workflow="$(printf '%s\n' "$entry" | sed -E 's/.* workflow=([^ ]+) sha=.*/\1/')"
    entry_sha="$(printf '%s\n' "$entry" | sed -E 's/.* sha=([^ ]+) -->$/\1/')"
    completion_id="$(printf '%s\n' "$completion" | sed -E 's/.* id=([^ ]+) slot=.*/\1/')"
    completion_slot="$(printf '%s\n' "$completion" | sed -E 's/.* slot=([0-9]{3}) workflow=.*/\1/')"
    completion_workflow="$(printf '%s\n' "$completion" | sed -E 's/.* workflow=([^ ]+) outcome=.*/\1/')"
    completion_outcome="$(printf '%s\n' "$completion" | sed -E 's/.* outcome=([^ ]+) source-commit=.*/\1/')"
    source_review_commit="$(printf '%s\n' "$completion" | sed -E 's/.* source-commit=([^ ]+) split-commit=.*/\1/')"
    split_review_commit="$(printf '%s\n' "$completion" | sed -E 's/.* split-commit=([^ ]+) -->$/\1/')"
    expected_slot="$(printf '%03d' "$((index + 2))")"

    test "$entry_slot" = "$expected_slot"
    test "$completion_id" = "$entry_id"
    test "$completion_slot" = "$entry_slot"
    test "$completion_workflow" = "$entry_workflow"
    entry_ids+=("$entry_id")
    entry_slots+=("$entry_slot")
    entry_kinds+=("$entry_kind")
    entry_modes+=("$entry_mode")
    entry_workflows+=("$entry_workflow")

    case "$entry_kind" in
        foundation)
            foundation_count=$((foundation_count + 1))
            test "$entry_sha" = "-"
            test "$entry_mode" = "workflow"
            ;;
        commit)
            test "$entry_sha" != "-"
            manifest_shas+=("$entry_sha")
            ;;
        residual)
            residual_count=$((residual_count + 1))
            test "$entry_sha" = "-"
            test "$entry_mode" = "workflow"
            ;;
    esac

    case "$entry_mode" in
        workflow)
            actual_workflow_count=$((actual_workflow_count + 1))
            test "$entry_workflow" != "-"
            workflow="${entry_workflow#@/}"
            test -f "$workflow"
            filename="${workflow#workflows/}"
            test "${filename%%_*}" = "$entry_slot"
            grep -Fq "$split_implementation" "$workflow"
            workflow_artifacts+=("$entry_workflow")

            case "$entry_kind" in
                foundation)
                    test "$(grep -Fxc 'id = "review_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_split_adjustment"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_split_adjustment"' "$workflow")" -eq 1
                    ;;
                commit)
                    grep -Fq "$entry_sha" "$workflow"
                    abbreviation="${entry_sha:0:12}"
                    printf '%s\n' "$workflow" |
                        grep -Eq "^workflows/${entry_slot}_${abbreviation}_[a-z0-9_-]+\\.donna\\.md$"
                    test "$(grep -Fxc 'id = "review_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_split_adjustment"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_split_adjustment"' "$workflow")" -eq 1
                    ;;
                residual)
                    test "$(grep -Fxc 'id = "verify_residual_inventory"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_residual_reports"' "$workflow")" -eq 1
                    ;;
            esac

            case "$entry_kind:$completion_outcome" in
                foundation:source-and-split|foundation:source-only|foundation:split-only|foundation:no-change)
                    ;;
                commit:source-and-split|commit:source-only|commit:split-only|commit:no-change)
                    ;;
                residual:reports)
                    ;;
                *)
                    exit 1
                    ;;
            esac
            ;;
        reference-only)
            actual_reference_count=$((actual_reference_count + 1))
            test "$entry_kind" = "commit"
            test "$entry_workflow" = "-"
            test "$completion_outcome" = "reference-only"
            test "$source_review_commit" = "none"
            test "$split_review_commit" = "none"
            test "$(find workflows -maxdepth 1 -type f -name "${entry_slot}_*.donna.md" -print | wc -l)" -eq 0
            reference_slots+=("$entry_slot")
            reference_shas+=("$entry_sha")
            ;;
    esac

    case "$completion_outcome" in
        source-and-split)
            test "$source_review_commit" != "none"
            test "$split_review_commit" != "none"
            ;;
        source-only)
            test "$source_review_commit" != "none"
            test "$split_review_commit" = "none"
            ;;
        split-only)
            test "$source_review_commit" = "none"
            test "$split_review_commit" != "none"
            ;;
        no-change|reference-only)
            test "$source_review_commit" = "none"
            test "$split_review_commit" = "none"
            ;;
        reports)
            test "$source_review_commit" = "none"
            ;;
    esac

    if test "$source_review_commit" != "none"; then
        git cat-file -e "$source_review_commit^{commit}"
        git merge-base --is-ancestor "$source_review_commit" HEAD
        review_commits+=("$source_review_commit")
        review_commit_count=$((review_commit_count + 1))
    fi
    if test "$split_review_commit" != "none"; then
        git cat-file -e "$split_review_commit^{commit}"
        git merge-base --is-ancestor "$split_review_commit" HEAD
        review_commits+=("$split_review_commit")
        review_commit_count=$((review_commit_count + 1))
    fi
    if test "$source_review_commit" != "none" && test "$split_review_commit" != "none"; then
        test "$source_review_commit" != "$split_review_commit"
        git merge-base --is-ancestor "$source_review_commit" "$split_review_commit"
    fi
done

test "$(
    printf '%s\n' "${review_commits[@]}" |
        sed '/^$/d' |
        sort -u |
        wc -l
)" -eq "$review_commit_count"
for ((index = 1; index < review_commit_count; index += 1)); do
    git merge-base --is-ancestor "${review_commits[$((index - 1))]}" "${review_commits[$index]}"
done

test "$actual_workflow_count" -eq "$workflow_count"
test "$actual_reference_count" -eq "$reference_count"
test "$foundation_count" -eq "$expected_foundations"
test "$residual_count" -eq 1
test "$(printf '%s\n' "${entry_ids[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${entry_slots[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${workflow_artifacts[@]}" | sort -u | wc -l)" -eq "$workflow_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | sort -u | wc -l)" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | cut -c1-12 | sort -u | wc -l)" -eq "$commit_count"

if test "$scope" = "post-sync"; then
    test "${entry_kinds[0]}" = "foundation"
    test "${entry_modes[0]}" = "workflow"
    test "${entry_workflows[0]}" = "@/workflows/002_foundation_core_api.donna.md"
    grep -Fq "$sync_merge" "${entry_workflows[0]#@/}"
    grep -Fq "$base" "${entry_workflows[0]#@/}"
    for ((index = 1; index <= commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
    done
    test "${entry_kinds[$((entry_count - 1))]}" = "residual"
    test "${reference_slots[*]}" = "005 018 024 025 026 029 030"
    test "${reference_shas[*]}" = "65f70a46c4cf44e2717744aeafbf8acbe83b0378 b8a262b4256966d60e1e8452ebde7a1bf471b4af 4a06d99949114b2804b8d34492a288872fb611ed 6601232c69b44c507f2e9f63d256836e655f7973 96d38eb7449e4d38c0d6a3fffd66e17f145f669f 0b3c3b55403f2bb7342a37236dd30f3d9b3544e9 20e76c9b9b28c0ec98faa63c3f9382c1347301b9"
    test "$(grep -Fxc '<!-- donna-entry kind=commit id=post-sync-015-4540aca9-generated-missing-docs slot=017 mode=workflow workflow=@/workflows/017_4540aca9182c_generated_missing_docs.donna.md sha=4540aca9182c6708fe9bda0b8fc33d2cf8c13e57 -->' "$manifest")" -eq 1
else
    for ((index = 0; index < commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
        test "${entry_modes[$index]}" = "workflow"
    done
    test "${entry_kinds[$commit_count]}" = "residual"
fi

mapfile -t expected_shas < <(git -C "$source_repo" rev-list --reverse --topo-order "$source_range")
test "$(printf '%s\n' "${expected_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
for index in "${!expected_shas[@]}"; do
    test "${manifest_shas[$index]}" = "${expected_shas[$index]}"
done

mapfile -t numbered_files < <(find workflows -maxdepth 1 -type f -name '[0-9][0-9][0-9]_*.donna.md' -print | sort)
numbered_count="$(printf '%s\n' "${numbered_files[@]}" | wc -l)"
test "$numbered_count" -eq "$((workflow_count + 1))"
test "${numbered_files[0]}" = "workflows/001_bec-ecoli-core-migration.donna.md"
expected_numbered_files=("workflows/001_bec-ecoli-core-migration.donna.md")
for artifact in "${workflow_artifacts[@]}"; do
    expected_numbered_files+=("${artifact#@/}")
done
for index in "${!numbered_files[@]}"; do
    test "${numbered_files[$index]}" = "${expected_numbered_files[$index]}"
done

split_marker_pattern='^<!-- split-path status=(A|D|M) path=[^ ]+ class=(required|recommended|equivalent|excluded|decision) -->$'
mapfile -t split_markers < <(grep -E "$split_marker_pattern" "$residual_report")
test "$(grep -Ec "$split_marker_pattern" "$residual_report")" -eq 99
test "$(grep -c '^<!-- split-path ' "$residual_report")" -eq 99
split_inventory_hash="$(
    printf '%s\n' "${split_markers[@]}" |
        sed -E 's/^<!-- split-path status=([ADM]) path=([^ ]+) class=[^ ]+ -->$/\1\t\2/' |
        sha256sum |
        awk '{print $1}'
)"
test "$split_inventory_hash" = "28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38"

for required_item in \
    "R/identifiers.R" \
    "R/matrix_adapter.R" \
    "R/registry.R" \
    "R/step_result.R" \
    "tests/testthat/helper-source-root.R" \
    "tests/testthat/test-identifiers.R" \
    "tests/testthat/test-lineage.R" \
    "tests/testthat/test-matrix_adapter.R" \
    "tests/testthat/test-registry.R" \
    "tests/testthat/test-step_result.R" \
    "tests/testthat/test-symbol-ownership.R" \
    "pb_apply_matrix_method" \
    "pb_step_result" \
    "pb_unregister_steps"; do
    grep -Fq "$required_item" "$residual_report"
done

bench_apis=(
    "correct_with_BERT"
    "correct_with_NormAE"
    "correct_with_PLSDA_batch"
    "correct_with_RUVIII_C"
    "correct_with_mComBat"
    "correct_with_omicsGMF"
    "estimate_omicsGMF_rank"
    "imputeMissForest"
    "imputeMissForest.ProBatchFeatures"
    "imputeMissForest_df"
    "imputeMissForest_dm"
    "missForestImpute"
    "imputePRONE"
    "imputePRONE_df"
    "imputePRONE_dm"
    "impute_and_correct_with_omicsGMF"
    "impute_with_omicsGMF"
    "impute_with_omicsGMF.ProBatchFeatures"
    "calculate_classification_metrics"
    "calculate_variance_partition"
    "prepare_variance_partition_df"
    "plot_variance_partition"
    "plot_variance_partition.df"
    "plot_intragroup_variation"
    "plot_TSNE"
    "plot_UMAP"
)
for bench_api in "${bench_apis[@]}"; do
    grep -Fxc -- "<!-- bench-api symbol=${bench_api} disposition=companion -->" "$residual_report"
done
test "$(grep -c '^<!-- bench-api ' "$residual_report")" -eq 26

actionable_keys="$(
    printf '%s\n' "${split_markers[@]}" |
        sed -nE 's/^<!-- split-path status=[ADM] path=([^ ]+) class=(required|recommended|decision) -->$/\2\t\1/p' |
        sort
)"
actionable_count="$(printf '%s\n' "$actionable_keys" | sed '/^$/d' | wc -l)"
plan_marker_pattern='^<!-- plan-path path=[^ ]+ class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$'
mapfile -t plan_markers < <(grep -E "$plan_marker_pattern" "$remaining_plan" || true)

if test "$actionable_count" -eq 0; then
    test "$(grep -Fxc '<!-- plan-empty -->' "$remaining_plan")" -eq 1
    test "$(grep -c '^<!-- plan-path ' "$remaining_plan" || true)" -eq 0
else
    test "$(grep -Fxc '<!-- plan-empty -->' "$remaining_plan" || true)" -eq 0
    test "$(grep -Ec "$plan_marker_pattern" "$remaining_plan")" -eq "$actionable_count"
    test "$(grep -c '^<!-- plan-path ' "$remaining_plan")" -eq "$actionable_count"

    plan_keys="$(
        printf '%s\n' "${plan_markers[@]}" |
            sed -E 's/^<!-- plan-path path=([^ ]+) class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$/\2\t\1/' |
            sort
    )"
    test "$plan_keys" = "$actionable_keys"

    mapfile -t plan_orders < <(
        printf '%s\n' "${plan_markers[@]}" |
            sed -E 's/.* order=([1-9][0-9]*) reason=.*/\1/' |
            sort -n
    )
    for index in "${!plan_orders[@]}"; do
        test "${plan_orders[$index]}" -eq "$((index + 1))"
    done
fi

donna -p llm validate --all
depmesh -p llm relations
depmesh -p llm dependencies --relation governs @/specs/general/workflows.md
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finalize Workflow Catalog

```toml donna
id = "finalize_catalog"
kind = "donna.lib.request_action"
```

1. Before editing the governed catalog, query all Depmesh relations for `{{ donna.lib.path("@/specs/general/workflows.md") }}` and every numbered workflow, and inspect every returned dependency.
2. Use the manifest and validated progress records to update the parent, every completed executable foundation/commit child, and the residual child from `implementation in progress` to `implemented` in the catalog. Reference-only records have no catalog artifacts or states. Do not change unrelated planned or implemented workflows.
3. Re-query Depmesh for the catalog and every numbered workflow, validate and render all workflows, and run generated-file-excluding diff checks.
4. Report only the final catalog/governance changes to the developer and give them time to inspect and manually commit this completion separately. The residual reports and progress records already passed `verify_completion`; do not edit them during this pause. Route any requested report, manifest, or progress correction through `repair_completion` so all completeness checks and this final review pause repeat.
5. Do not create or amend that commit, and do not complete this action request until the developer explicitly resumes.
6. If the developer accepts the final catalog and resumes, `{{ donna.lib.goto("verify_final_catalog") }}`.
7. If catalog or governance corrections are requested, `{{ donna.lib.goto("repair_final_catalog") }}`.
8. If report, manifest, or progress corrections are requested, `{{ donna.lib.goto("repair_completion") }}`.
9. If the completion state cannot be accepted, `{{ donna.lib.goto("blocked") }}`.

## Verify Final Catalog

```toml donna
id = "verify_final_catalog"
kind = "donna.lib.run_script"
save_stdout_to = "final_catalog_stdout"
save_stderr_to = "final_catalog_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_final_catalog"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

manifest="workflows/bec-ecoli-core-migration-manifest.md"
progress=".session/donna/bec-ecoli-core-migration-progress.md"
residual_report="workflows/reports/bec-ecoli-core-residual-review.md"
remaining_plan="workflows/reports/bec-ecoli-core-remaining-change-plan.md"
catalog="specs/general/workflows.md"
source_repo="/home/yuliya/repos/other/proBatch"
base="ba6ee246eace090e71baa7aba302ca64e76ddb32"
sync_merge="5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab"
tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"
split_implementation="29a7478dc7deea846a2c1ff1abd25a881e6f87db"

test -f "$manifest"
test -f "$progress"
test -s "$residual_report"
test -s "$remaining_plan"

header_pattern='^<!-- donna-migration scope=(post-sync|full-difference) last-prefix=[0-9]{3} entries=[0-9]+ workflows=[0-9]+ commits=[0-9]+ references=[0-9]+ -->$'
mapfile -t headers < <(grep -E "$header_pattern" "$manifest")
test "$(grep -Ec "$header_pattern" "$manifest")" -eq 1
read -r scope last_prefix entry_count workflow_count commit_count reference_count < <(
    printf '%s\n' "${headers[0]}" |
        sed -E 's/^<!-- donna-migration scope=([^ ]+) last-prefix=([0-9]{3}) entries=([0-9]+) workflows=([0-9]+) commits=([0-9]+) references=([0-9]+) -->$/\1 \2 \3 \4 \5 \6/'
)

case "$scope" in
    post-sync)
        test "$last_prefix" = "032"
        test "$entry_count" -eq 31
        test "$workflow_count" -eq 24
        test "$commit_count" -eq 29
        test "$reference_count" -eq 7
        source_range="$sync_merge..$tip"
        expected_foundations=1
        ;;
    full-difference)
        test "$last_prefix" = "314"
        test "$entry_count" -eq 313
        test "$workflow_count" -eq 313
        test "$commit_count" -eq 312
        test "$reference_count" -eq 0
        source_range="$base..$tip"
        expected_foundations=0
        ;;
    *)
        exit 1
        ;;
esac

entry_pattern='^<!-- donna-entry kind=(foundation|commit|residual) id=[a-z0-9-]+ slot=[0-9]{3} mode=(workflow|reference-only) workflow=(-|@/workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md) sha=(-|[0-9a-f]{40}) -->$'
completion_pattern='^<!-- donna-complete id=[a-z0-9-]+ slot=[0-9]{3} workflow=(-|@/workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md) outcome=(source-and-split|source-only|split-only|no-change|reference-only|reports) source-commit=(none|[0-9a-f]{40}) split-commit=(none|[0-9a-f]{40}) -->$'
mapfile -t entries < <(grep -E "$entry_pattern" "$manifest")
mapfile -t completions < <(grep -E "$completion_pattern" "$progress")
test "$(grep -Ec "$entry_pattern" "$manifest")" -eq "$entry_count"
test "$(grep -c '^<!-- donna-entry ' "$manifest")" -eq "$entry_count"
test "$(grep -Ec "$completion_pattern" "$progress")" -eq "$entry_count"
test "$(grep -c '^<!-- donna-complete ' "$progress")" -eq "$entry_count"

entry_ids=()
entry_slots=()
entry_kinds=()
entry_modes=()
entry_workflows=()
workflow_artifacts=()
manifest_shas=()
reference_slots=()
reference_shas=()
foundation_count=0
residual_count=0
actual_workflow_count=0
actual_reference_count=0
review_commits=()
review_commit_count=0

for index in "${!entries[@]}"; do
    entry="${entries[$index]}"
    completion="${completions[$index]}"
    entry_kind="$(printf '%s\n' "$entry" | sed -E 's/^<!-- donna-entry kind=([^ ]+) .*/\1/')"
    entry_id="$(printf '%s\n' "$entry" | sed -E 's/.* id=([^ ]+) slot=.*/\1/')"
    entry_slot="$(printf '%s\n' "$entry" | sed -E 's/.* slot=([0-9]{3}) mode=.*/\1/')"
    entry_mode="$(printf '%s\n' "$entry" | sed -E 's/.* mode=([^ ]+) workflow=.*/\1/')"
    entry_workflow="$(printf '%s\n' "$entry" | sed -E 's/.* workflow=([^ ]+) sha=.*/\1/')"
    entry_sha="$(printf '%s\n' "$entry" | sed -E 's/.* sha=([^ ]+) -->$/\1/')"
    completion_id="$(printf '%s\n' "$completion" | sed -E 's/.* id=([^ ]+) slot=.*/\1/')"
    completion_slot="$(printf '%s\n' "$completion" | sed -E 's/.* slot=([0-9]{3}) workflow=.*/\1/')"
    completion_workflow="$(printf '%s\n' "$completion" | sed -E 's/.* workflow=([^ ]+) outcome=.*/\1/')"
    completion_outcome="$(printf '%s\n' "$completion" | sed -E 's/.* outcome=([^ ]+) source-commit=.*/\1/')"
    source_review_commit="$(printf '%s\n' "$completion" | sed -E 's/.* source-commit=([^ ]+) split-commit=.*/\1/')"
    split_review_commit="$(printf '%s\n' "$completion" | sed -E 's/.* split-commit=([^ ]+) -->$/\1/')"
    expected_slot="$(printf '%03d' "$((index + 2))")"

    test "$entry_slot" = "$expected_slot"
    test "$completion_id" = "$entry_id"
    test "$completion_slot" = "$entry_slot"
    test "$completion_workflow" = "$entry_workflow"
    entry_ids+=("$entry_id")
    entry_slots+=("$entry_slot")
    entry_kinds+=("$entry_kind")
    entry_modes+=("$entry_mode")
    entry_workflows+=("$entry_workflow")

    case "$entry_kind" in
        foundation)
            foundation_count=$((foundation_count + 1))
            test "$entry_sha" = "-"
            test "$entry_mode" = "workflow"
            ;;
        commit)
            test "$entry_sha" != "-"
            manifest_shas+=("$entry_sha")
            ;;
        residual)
            residual_count=$((residual_count + 1))
            test "$entry_sha" = "-"
            test "$entry_mode" = "workflow"
            ;;
    esac

    case "$entry_mode" in
        workflow)
            actual_workflow_count=$((actual_workflow_count + 1))
            test "$entry_workflow" != "-"
            workflow="${entry_workflow#@/}"
            test -f "$workflow"
            filename="${workflow#workflows/}"
            test "${filename%%_*}" = "$entry_slot"
            grep -Fq "$split_implementation" "$workflow"
            workflow_artifacts+=("$entry_workflow")

            case "$entry_kind" in
                foundation)
                    test "$(grep -Fxc 'id = "review_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_split_adjustment"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_split_adjustment"' "$workflow")" -eq 1
                    ;;
                commit)
                    grep -Fq "$entry_sha" "$workflow"
                    abbreviation="${entry_sha:0:12}"
                    printf '%s\n' "$workflow" |
                        grep -Eq "^workflows/${entry_slot}_${abbreviation}_[a-z0-9_-]+\\.donna\\.md$"
                    test "$(grep -Fxc 'id = "review_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_bec_port"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_split_adjustment"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "reverify_split_adjustment"' "$workflow")" -eq 1
                    ;;
                residual)
                    test "$(grep -Fxc 'id = "verify_residual_inventory"' "$workflow")" -eq 1
                    test "$(grep -Fxc 'id = "review_residual_reports"' "$workflow")" -eq 1
                    ;;
            esac

            case "$entry_kind:$completion_outcome" in
                foundation:source-and-split|foundation:source-only|foundation:split-only|foundation:no-change)
                    ;;
                commit:source-and-split|commit:source-only|commit:split-only|commit:no-change)
                    ;;
                residual:reports)
                    ;;
                *)
                    exit 1
                    ;;
            esac
            ;;
        reference-only)
            actual_reference_count=$((actual_reference_count + 1))
            test "$entry_kind" = "commit"
            test "$entry_workflow" = "-"
            test "$completion_outcome" = "reference-only"
            test "$source_review_commit" = "none"
            test "$split_review_commit" = "none"
            test "$(find workflows -maxdepth 1 -type f -name "${entry_slot}_*.donna.md" -print | wc -l)" -eq 0
            reference_slots+=("$entry_slot")
            reference_shas+=("$entry_sha")
            ;;
    esac

    case "$completion_outcome" in
        source-and-split)
            test "$source_review_commit" != "none"
            test "$split_review_commit" != "none"
            ;;
        source-only)
            test "$source_review_commit" != "none"
            test "$split_review_commit" = "none"
            ;;
        split-only)
            test "$source_review_commit" = "none"
            test "$split_review_commit" != "none"
            ;;
        no-change|reference-only)
            test "$source_review_commit" = "none"
            test "$split_review_commit" = "none"
            ;;
        reports)
            test "$source_review_commit" = "none"
            ;;
    esac

    if test "$source_review_commit" != "none"; then
        git cat-file -e "$source_review_commit^{commit}"
        git merge-base --is-ancestor "$source_review_commit" HEAD
        review_commits+=("$source_review_commit")
        review_commit_count=$((review_commit_count + 1))
    fi
    if test "$split_review_commit" != "none"; then
        git cat-file -e "$split_review_commit^{commit}"
        git merge-base --is-ancestor "$split_review_commit" HEAD
        review_commits+=("$split_review_commit")
        review_commit_count=$((review_commit_count + 1))
    fi
    if test "$source_review_commit" != "none" && test "$split_review_commit" != "none"; then
        test "$source_review_commit" != "$split_review_commit"
        git merge-base --is-ancestor "$source_review_commit" "$split_review_commit"
    fi
done

test "$(
    printf '%s\n' "${review_commits[@]}" |
        sed '/^$/d' |
        sort -u |
        wc -l
)" -eq "$review_commit_count"
for ((index = 1; index < review_commit_count; index += 1)); do
    git merge-base --is-ancestor "${review_commits[$((index - 1))]}" "${review_commits[$index]}"
done

test "$actual_workflow_count" -eq "$workflow_count"
test "$actual_reference_count" -eq "$reference_count"
test "$foundation_count" -eq "$expected_foundations"
test "$residual_count" -eq 1
test "$(printf '%s\n' "${entry_ids[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${entry_slots[@]}" | sort -u | wc -l)" -eq "$entry_count"
test "$(printf '%s\n' "${workflow_artifacts[@]}" | sort -u | wc -l)" -eq "$workflow_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | sort -u | wc -l)" -eq "$commit_count"
test "$(printf '%s\n' "${manifest_shas[@]}" | cut -c1-12 | sort -u | wc -l)" -eq "$commit_count"

if test "$scope" = "post-sync"; then
    test "${entry_kinds[0]}" = "foundation"
    test "${entry_modes[0]}" = "workflow"
    test "${entry_workflows[0]}" = "@/workflows/002_foundation_core_api.donna.md"
    grep -Fq "$sync_merge" "${entry_workflows[0]#@/}"
    grep -Fq "$base" "${entry_workflows[0]#@/}"
    for ((index = 1; index <= commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
    done
    test "${entry_kinds[$((entry_count - 1))]}" = "residual"
    test "${reference_slots[*]}" = "005 018 024 025 026 029 030"
    test "${reference_shas[*]}" = "65f70a46c4cf44e2717744aeafbf8acbe83b0378 b8a262b4256966d60e1e8452ebde7a1bf471b4af 4a06d99949114b2804b8d34492a288872fb611ed 6601232c69b44c507f2e9f63d256836e655f7973 96d38eb7449e4d38c0d6a3fffd66e17f145f669f 0b3c3b55403f2bb7342a37236dd30f3d9b3544e9 20e76c9b9b28c0ec98faa63c3f9382c1347301b9"
    test "$(grep -Fxc '<!-- donna-entry kind=commit id=post-sync-015-4540aca9-generated-missing-docs slot=017 mode=workflow workflow=@/workflows/017_4540aca9182c_generated_missing_docs.donna.md sha=4540aca9182c6708fe9bda0b8fc33d2cf8c13e57 -->' "$manifest")" -eq 1
else
    for ((index = 0; index < commit_count; index += 1)); do
        test "${entry_kinds[$index]}" = "commit"
        test "${entry_modes[$index]}" = "workflow"
    done
    test "${entry_kinds[$commit_count]}" = "residual"
fi

mapfile -t expected_shas < <(git -C "$source_repo" rev-list --reverse --topo-order "$source_range")
test "$(printf '%s\n' "${expected_shas[@]}" | grep -Ec '^[0-9a-f]{40}$')" -eq "$commit_count"
for index in "${!expected_shas[@]}"; do
    test "${manifest_shas[$index]}" = "${expected_shas[$index]}"
done

split_marker_pattern='^<!-- split-path status=(A|D|M) path=[^ ]+ class=(required|recommended|equivalent|excluded|decision) -->$'
mapfile -t split_markers < <(grep -E "$split_marker_pattern" "$residual_report")
test "$(grep -Ec "$split_marker_pattern" "$residual_report")" -eq 99
test "$(grep -c '^<!-- split-path ' "$residual_report")" -eq 99
split_inventory_hash="$(
    printf '%s\n' "${split_markers[@]}" |
        sed -E 's/^<!-- split-path status=([ADM]) path=([^ ]+) class=[^ ]+ -->$/\1\t\2/' |
        sha256sum |
        awk '{print $1}'
)"
test "$split_inventory_hash" = "28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38"

for required_item in \
    "R/identifiers.R" \
    "R/matrix_adapter.R" \
    "R/registry.R" \
    "R/step_result.R" \
    "tests/testthat/helper-source-root.R" \
    "tests/testthat/test-identifiers.R" \
    "tests/testthat/test-lineage.R" \
    "tests/testthat/test-matrix_adapter.R" \
    "tests/testthat/test-registry.R" \
    "tests/testthat/test-step_result.R" \
    "tests/testthat/test-symbol-ownership.R" \
    "pb_apply_matrix_method" \
    "pb_step_result" \
    "pb_unregister_steps"; do
    grep -Fq "$required_item" "$residual_report"
done

bench_apis=(
    "correct_with_BERT"
    "correct_with_NormAE"
    "correct_with_PLSDA_batch"
    "correct_with_RUVIII_C"
    "correct_with_mComBat"
    "correct_with_omicsGMF"
    "estimate_omicsGMF_rank"
    "imputeMissForest"
    "imputeMissForest.ProBatchFeatures"
    "imputeMissForest_df"
    "imputeMissForest_dm"
    "missForestImpute"
    "imputePRONE"
    "imputePRONE_df"
    "imputePRONE_dm"
    "impute_and_correct_with_omicsGMF"
    "impute_with_omicsGMF"
    "impute_with_omicsGMF.ProBatchFeatures"
    "calculate_classification_metrics"
    "calculate_variance_partition"
    "prepare_variance_partition_df"
    "plot_variance_partition"
    "plot_variance_partition.df"
    "plot_intragroup_variation"
    "plot_TSNE"
    "plot_UMAP"
)
for bench_api in "${bench_apis[@]}"; do
    grep -Fxc -- "<!-- bench-api symbol=${bench_api} disposition=companion -->" "$residual_report"
done
test "$(grep -c '^<!-- bench-api ' "$residual_report")" -eq 26

actionable_keys="$(
    printf '%s\n' "${split_markers[@]}" |
        sed -nE 's/^<!-- split-path status=[ADM] path=([^ ]+) class=(required|recommended|decision) -->$/\2\t\1/p' |
        sort
)"
actionable_count="$(printf '%s\n' "$actionable_keys" | sed '/^$/d' | wc -l)"
plan_marker_pattern='^<!-- plan-path path=[^ ]+ class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$'
mapfile -t plan_markers < <(grep -E "$plan_marker_pattern" "$remaining_plan" || true)

if test "$actionable_count" -eq 0; then
    test "$(grep -Fxc '<!-- plan-empty -->' "$remaining_plan")" -eq 1
    test "$(grep -c '^<!-- plan-path ' "$remaining_plan" || true)" -eq 0
else
    test "$(grep -Fxc '<!-- plan-empty -->' "$remaining_plan" || true)" -eq 0
    test "$(grep -Ec "$plan_marker_pattern" "$remaining_plan")" -eq "$actionable_count"
    test "$(grep -c '^<!-- plan-path ' "$remaining_plan")" -eq "$actionable_count"

    plan_keys="$(
        printf '%s\n' "${plan_markers[@]}" |
            sed -E 's/^<!-- plan-path path=([^ ]+) class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$/\2\t\1/' |
            sort
    )"
    test "$plan_keys" = "$actionable_keys"

    mapfile -t plan_orders < <(
        printf '%s\n' "${plan_markers[@]}" |
            sed -E 's/.* order=([1-9][0-9]*) reason=.*/\1/' |
            sort -n
    )
    for index in "${!plan_orders[@]}"; do
        test "${plan_orders[$index]}" -eq "$((index + 1))"
    done
fi

mapfile -t numbered_files < <(find workflows -maxdepth 1 -type f -name '[0-9][0-9][0-9]_*.donna.md' -print | sort)
numbered_count="$(printf '%s\n' "${numbered_files[@]}" | wc -l)"
test "$numbered_count" -eq "$((workflow_count + 1))"
test "${numbered_files[0]}" = "workflows/001_bec-ecoli-core-migration.donna.md"
test "$(printf '%s\n' "${numbered_files[@]}" | grep -Ec '^workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md$')" -eq "$numbered_count"
expected_numbered_files=("workflows/001_bec-ecoli-core-migration.donna.md")
for artifact in "${workflow_artifacts[@]}"; do
    expected_numbered_files+=("${artifact#@/}")
done
for index in "${!numbered_files[@]}"; do
    test "${numbered_files[$index]}" = "${expected_numbered_files[$index]}"
done

catalog_pattern='^- Artifact: `\./workflows/[0-9]{3}_[a-z0-9_-]+\.donna\.md`$'
test "$(grep -Ec "$catalog_pattern" "$catalog")" -eq "$numbered_count"
governs_output="$(depmesh -p llm dependencies --relation governs @/specs/general/workflows.md)"
for workflow in "${numbered_files[@]}"; do
    artifact="@/${workflow}"
    catalog_line="- Artifact: \`./${workflow}\`"
    dependency_output="$(depmesh -p llm dependencies "$artifact")"

    test "$(grep -Fxc -- "$catalog_line" "$catalog")" -eq 1
    awk -v target="$catalog_line" '
        $0 == target {
            getline
            found = ($0 == "- State: implemented")
        }
        END {
            exit(found ? 0 : 1)
        }
    ' "$catalog"
    printf '%s\n' "$dependency_output" | grep -Fxq -- "- @/specs/behavior/files_relations.md"
    printf '%s\n' "$dependency_output" | grep -Fxq -- "- @/specs/general/workflows.md"
    test "$(printf '%s\n' "$dependency_output" | grep -c '^- @/specs/')" -eq 2
    printf '%s\n' "$governs_output" | grep -Fxq -- "- $artifact"

    donna -p llm validate "$artifact"
    donna -p llm render "$artifact" --mode view >/dev/null
done

donna -p llm validate --all
depmesh -p llm relations
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair Final Catalog

```toml donna
id = "repair_final_catalog"
kind = "donna.lib.request_action"
```

The final catalog or post-review completeness gate failed, or the developer requested a catalog/governance correction.

Final-gate standard output:

```text
{{ donna.lib.task_variable("final_catalog_stdout") }}
```

Final-gate standard error:

```text
{{ donna.lib.task_variable("final_catalog_stderr") }}
```

1. Diagnose the exact catalog, governance, workflow, progress, or residual-report failure. Before editing any governed artifact, query all its Depmesh relations and inspect the returned dependencies.
2. If a catalog, governance, or workflow artifact must change, make the narrow correction and `{{ donna.lib.goto("finalize_catalog") }}` so dependency verification and the developer's manual catalog review/commit pause repeat.
3. If a manifest, progress record, or residual report must change, make the narrow correction and `{{ donna.lib.goto("verify_completion") }}` so the complete residual gate and final catalog pause repeat.
4. Only when diagnosis confirms that no artifact changed and the failure was transient, `{{ donna.lib.goto("verify_final_catalog") }}`.
5. If correction requires changing developer-reviewed package behavior or an unresolved ownership decision, `{{ donna.lib.goto("blocked") }}`.

## Repair Completion

```toml donna
id = "repair_completion"
kind = "donna.lib.request_action"
```

The final completeness gate failed.

Standard output:

```text
{{ donna.lib.task_variable("completion_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("completion_stderr") }}
```

1. Diagnose missing workflow records, reports, or residual validation failures without changing developer-reviewed package behavior. Before editing a governed artifact, query all its Depmesh relations and inspect the returned dependencies.
2. After any correction—or after restoring a transient input without an artifact edit—`{{ donna.lib.goto("verify_completion") }}`. A successful recheck must repeat the final catalog review/commit pause before the workflow can finish.
3. If developer input or a package-behavior change is required, `{{ donna.lib.goto("blocked") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The selected BEC history has been reviewed in order, every core-relevant source and stopped-split adjustment has received its own maintainer-controlled review opportunity, and the remaining split recommendations are documented. Report the selected history scope, completed numbered workflows, developer commits or no-change decisions, verification results, residual report paths, and unresolved API decisions.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The migration cannot continue safely. Report the completed numbered workflows, exact source SHA or workflow at the stop point, the unresolved history/ownership/API decision, preserved working-tree state, and the explicit developer input required.
