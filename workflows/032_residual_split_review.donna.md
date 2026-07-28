# Review Residual Stopped-Split Evidence

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_residual_evidence"
```

Audit every non-generated path in the rejected stopped split after the ordered BEC migration, produce a residual review and ordered change plan, and pause for maintainer review.

## Inspect Residual Evidence

```toml donna
id = "inspect_residual_evidence"
kind = "donna.lib.request_action"
```

1. Read the `residual-split-review` entry in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}` and the completed migration progress in `{{ donna.lib.path("@/.session/donna/bec-ecoli-core-migration-progress.md") }}`.
2. Read `/home/yuliya/repos/other/proBatch/archive/stopped-probatch-split-2026-07-27/API_FUNCTIONS.md` completely.
3. Read the stopped provenance note only with `git -C /home/yuliya/repos/other/proBatch-core-split show a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md`.
4. Account for the non-generated baseline-to-split change `ba6ee246eace090e71baa7aba302ca64e76ddb32..29a7478dc7deea846a2c1ff1abd25a881e6f87db`, including mode-only preflight `49cee7cc978fbb149c262a5a783face32dd1d135`.
5. Build the exhaustive path inventory with the pinned Git-object diff `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92..29a7478dc7deea846a2c1ff1abd25a881e6f87db`. Use only `git -C ... show`, `diff`, `ls-tree`, or `grep` against full IDs.
6. Never inspect either external working tree, write to an external repository, or run checkout, switch, reset, stash, or clean. Never read, enumerate, compare, or validate `man/**` or generated `NAMESPACE`; exclude them in every Git command.
7. Explicitly assess:
   - split-only `R/identifiers.R`, `R/matrix_adapter.R`, `R/registry.R`, and `R/step_result.R`;
   - `tests/testthat/helper-source-root.R`, `test-identifiers.R`, `test-lineage.R`, `test-matrix_adapter.R`, `test-registry.R`, `test-step_result.R`, and `test-symbol-ownership.R`;
   - new APIs `pb_apply_matrix_method`, `pb_step_result`, and `pb_unregister_steps`;
   - all 26 companion-owned exports, provider removals, duplicate-definition consolidation, and retained/deleted tests.
8. If the evidence is complete, `{{ donna.lib.goto("draft_residual_reports") }}`.
9. If a pinned identity or ownership boundary is unresolved, `{{ donna.lib.goto("blocked") }}`.

## Draft Residual Reports

```toml donna
id = "draft_residual_reports"
kind = "donna.lib.request_action"
```

1. Query all Depmesh relations for both report artifacts before editing.
2. Create or update:
   - `{{ donna.lib.path("@/workflows/reports/bec-ecoli-core-residual-review.md") }}`
   - `{{ donna.lib.path("@/workflows/reports/bec-ecoli-core-remaining-change-plan.md") }}`
3. Preserve unrelated changes and use `apply_patch`. Do not stage, commit, amend, or draft a commit message.
4. In the residual review, account for all 99 immutable paths exactly once and in Git name-status order with `<!-- split-path status=<A|D|M> path=<path> class=<required|recommended|equivalent|excluded|decision> -->`.
5. Explain each classification with the source symbol/path, destination, originating evidence, reason, dependency effect, focused tests, ordering, and explicit exclusion where applicable.
6. Add exactly one `<!-- bench-api symbol=<exact-name> disposition=companion -->` marker for each of the 26 companion exports listed in the ownership map.
7. In the remaining-change plan, add exactly one `<!-- plan-path path=<path> class=<required|recommended|decision> destination=<nonspace-target> order=<positive-integer> reason=<nonspace-reason-key> -->` marker for every actionable residual path and explain it. If there are no actionable paths, use exactly `<!-- plan-empty -->`.
8. Record audit commands and results in `.session/donna/032-residual-split-review-notes.md`.
9. When both reports are complete, `{{ donna.lib.goto("verify_residual_inventory") }}`.
10. If a classification needs maintainer input before a complete inventory can be written, use class `decision`; if even that is unsafe, `{{ donna.lib.goto("blocked") }}`.

## Verify Residual Inventory

```toml donna
id = "verify_residual_inventory"
kind = "donna.lib.run_script"
save_stdout_to = "residual_verify_stdout"
save_stderr_to = "residual_verify_stderr"
goto_on_success = "review_residual_reports"
goto_on_failure = "repair_residual_reports"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

review="workflows/reports/bec-ecoli-core-residual-review.md"
plan="workflows/reports/bec-ecoli-core-remaining-change-plan.md"
split_repo="/home/yuliya/repos/other/proBatch-core-split"
source_tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"
split_commit="29a7478dc7deea846a2c1ff1abd25a881e6f87db"
expected_checksum="28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38"

test -f "$review"
test -f "$plan"

expected="$(mktemp)"
actual="$(mktemp)"
trap 'rm -f "$expected" "$actual"' EXIT

git -C "$split_repo" diff --name-status "$source_tip" "$split_commit" -- . \
    ':(exclude)man/**' ':(exclude)NAMESPACE' >"$expected"

test "$(wc -l <"$expected")" -eq 99
test "$(awk '$1 == "A" { count++ } END { print count + 0 }' "$expected")" -eq 11
test "$(awk '$1 == "D" { count++ } END { print count + 0 }' "$expected")" -eq 38
test "$(awk '$1 == "M" { count++ } END { print count + 0 }' "$expected")" -eq 50
test "$(sha256sum "$expected" | awk '{ print $1 }')" = "$expected_checksum"

marker_pattern='^<!-- split-path status=(A|D|M) path=[^ ]+ class=(required|recommended|equivalent|excluded|decision) -->$'
test "$(grep -Ec "$marker_pattern" "$review")" -eq 99
test "$(grep -c '^<!-- split-path ' "$review")" -eq 99
sed -E 's/^<!-- split-path status=([ADM]) path=([^ ]+) class=(required|recommended|equivalent|excluded|decision) -->$/\1\t\2/' \
    < <(grep -E "$marker_pattern" "$review") >"$actual"
cmp "$expected" "$actual"

bench_symbols=(
    correct_with_BERT
    correct_with_NormAE
    correct_with_PLSDA_batch
    correct_with_RUVIII_C
    correct_with_mComBat
    correct_with_omicsGMF
    estimate_omicsGMF_rank
    imputeMissForest
    imputeMissForest.ProBatchFeatures
    imputeMissForest_df
    imputeMissForest_dm
    missForestImpute
    imputePRONE
    imputePRONE_df
    imputePRONE_dm
    impute_and_correct_with_omicsGMF
    impute_with_omicsGMF
    impute_with_omicsGMF.ProBatchFeatures
    calculate_classification_metrics
    calculate_variance_partition
    prepare_variance_partition_df
    plot_variance_partition
    plot_variance_partition.df
    plot_intragroup_variation
    plot_TSNE
    plot_UMAP
)

test "$(grep -Ec '^<!-- bench-api symbol=[^ ]+ disposition=companion -->$' "$review")" -eq 26
test "$(grep -c '^<!-- bench-api ' "$review")" -eq 26
for symbol in "${bench_symbols[@]}"; do
    test "$(grep -Fxc "<!-- bench-api symbol=$symbol disposition=companion -->" "$review")" -eq 1
done

mapfile -t actionable < <(
    sed -n -E \
        's/^<!-- split-path status=[ADM] path=([^ ]+) class=(required|recommended|decision) -->$/\1\t\2/p' \
        "$review"
)
actionable_count="$(printf '%s\n' "${actionable[@]}" | sed '/^$/d' | wc -l)"
plan_pattern='^<!-- plan-path path=[^ ]+ class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$'
plan_count="$(grep -Ec "$plan_pattern" "$plan" || true)"
test "$(grep -c '^<!-- plan-path ' "$plan" || true)" -eq "$plan_count"

if test "$actionable_count" -eq 0; then
    test "$plan_count" -eq 0
    test "$(grep -Fxc '<!-- plan-empty -->' "$plan")" -eq 1
else
    test "$plan_count" -eq "$actionable_count"
    test "$(grep -Fxc '<!-- plan-empty -->' "$plan" || true)" -eq 0

    declare -A plan_classes=()
    mapfile -t plan_markers < <(grep -E "$plan_pattern" "$plan")
    for marker in "${plan_markers[@]}"; do
        path="$(printf '%s\n' "$marker" | sed -E 's/^<!-- plan-path path=([^ ]+) .*/\1/')"
        class="$(printf '%s\n' "$marker" | sed -E 's/.* class=([^ ]+) destination=.*/\1/')"
        test -z "${plan_classes[$path]+x}"
        plan_classes["$path"]="$class"
    done

    for item in "${actionable[@]}"; do
        path="${item%%$'\t'*}"
        class="${item#*$'\t'}"
        test "${plan_classes[$path]+x}" = x
        test "${plan_classes[$path]}" = "$class"
    done
fi

git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair Residual Reports

```toml donna
id = "repair_residual_reports"
kind = "donna.lib.request_action"
```

Residual inventory verification failed.

Standard output:

```text
{{ donna.lib.task_variable("residual_verify_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("residual_verify_stderr") }}
```

1. Repair only the two reports or their session evidence. Query Depmesh before editing either governed report.
2. Preserve all completed package behavior and the immutable 99-path order.
3. When repaired, `{{ donna.lib.goto("verify_residual_inventory") }}`.
4. If repair requires changing maintainer-reviewed package behavior, `{{ donna.lib.goto("blocked") }}`.

## Review Residual Reports

```toml donna
id = "review_residual_reports"
kind = "donna.lib.request_action"
```

This is a maintainer-controlled report review and commit pause.

1. Report the 99-path classification totals, all `decision` items, the ordered actionable plan, explicit exclusions, and verification result.
2. Give the maintainer time to inspect, edit, stage, and manually commit the reports if desired.
3. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
4. Do not complete this action request until the maintainer explicitly asks to resume and confirms either the report commit SHA or an explicit no-new-commit decision.
5. On explicit resume, record exactly one `<!-- residual-review commit=<40-hex-or-none> -->` marker in `.session/donna/032-residual-split-review-notes.md`.
6. If accepted, `{{ donna.lib.goto("reverify_residual_inventory") }}`.
7. If corrections are requested, `{{ donna.lib.goto("repair_residual_reports") }}`.
8. If rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

## Reverify Residual Inventory

```toml donna
id = "reverify_residual_inventory"
kind = "donna.lib.run_script"
save_stdout_to = "residual_reverify_stdout"
save_stderr_to = "residual_reverify_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_residual_reports"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

notes=".session/donna/032-residual-split-review-notes.md"
review="workflows/reports/bec-ecoli-core-residual-review.md"
plan="workflows/reports/bec-ecoli-core-remaining-change-plan.md"
split_repo="/home/yuliya/repos/other/proBatch-core-split"
source_tip="e2bb18547c73f1c471fc1afcb3facbd8bea5fa92"
split_commit="29a7478dc7deea846a2c1ff1abd25a881e6f87db"

test -f "$notes"
test "$(grep -Ec '^<!-- residual-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
test -f "$review"
test -f "$plan"

expected="$(mktemp)"
actual="$(mktemp)"
trap 'rm -f "$expected" "$actual"' EXIT
git -C "$split_repo" diff --name-status "$source_tip" "$split_commit" -- . \
    ':(exclude)man/**' ':(exclude)NAMESPACE' >"$expected"
test "$(wc -l <"$expected")" -eq 99
test "$(sha256sum "$expected" | awk '{ print $1 }')" = "28828d60c51178f042ca3f2389255bb69527ef18cc5aa4c5cdd8e4b687274b38"

marker_pattern='^<!-- split-path status=(A|D|M) path=[^ ]+ class=(required|recommended|equivalent|excluded|decision) -->$'
test "$(grep -Ec "$marker_pattern" "$review")" -eq 99
test "$(grep -c '^<!-- split-path ' "$review")" -eq 99
sed -E 's/^<!-- split-path status=([ADM]) path=([^ ]+) class=(required|recommended|equivalent|excluded|decision) -->$/\1\t\2/' \
    < <(grep -E "$marker_pattern" "$review") >"$actual"
cmp "$expected" "$actual"
bench_symbols=(
    correct_with_BERT
    correct_with_NormAE
    correct_with_PLSDA_batch
    correct_with_RUVIII_C
    correct_with_mComBat
    correct_with_omicsGMF
    estimate_omicsGMF_rank
    imputeMissForest
    imputeMissForest.ProBatchFeatures
    imputeMissForest_df
    imputeMissForest_dm
    missForestImpute
    imputePRONE
    imputePRONE_df
    imputePRONE_dm
    impute_and_correct_with_omicsGMF
    impute_with_omicsGMF
    impute_with_omicsGMF.ProBatchFeatures
    calculate_classification_metrics
    calculate_variance_partition
    prepare_variance_partition_df
    plot_variance_partition
    plot_variance_partition.df
    plot_intragroup_variation
    plot_TSNE
    plot_UMAP
)
test "$(grep -Ec '^<!-- bench-api symbol=[^ ]+ disposition=companion -->$' "$review")" -eq 26
test "$(grep -c '^<!-- bench-api ' "$review")" -eq 26
for symbol in "${bench_symbols[@]}"; do
    test "$(grep -Fxc "<!-- bench-api symbol=$symbol disposition=companion -->" "$review")" -eq 1
done

mapfile -t actionable < <(
    sed -n -E \
        's/^<!-- split-path status=[ADM] path=([^ ]+) class=(required|recommended|decision) -->$/\1\t\2/p' \
        "$review"
)
actionable_count="$(printf '%s\n' "${actionable[@]}" | sed '/^$/d' | wc -l)"
plan_pattern='^<!-- plan-path path=[^ ]+ class=(required|recommended|decision) destination=[^ ]+ order=[1-9][0-9]* reason=[^ ]+ -->$'
plan_count="$(grep -Ec "$plan_pattern" "$plan" || true)"
test "$(grep -c '^<!-- plan-path ' "$plan" || true)" -eq "$plan_count"
if test "$actionable_count" -eq 0; then
    test "$plan_count" -eq 0
    test "$(grep -Fxc '<!-- plan-empty -->' "$plan")" -eq 1
else
    test "$plan_count" -eq "$actionable_count"
    test "$(grep -Fxc '<!-- plan-empty -->' "$plan" || true)" -eq 0

    declare -A plan_classes=()
    mapfile -t plan_markers < <(grep -E "$plan_pattern" "$plan")
    for marker in "${plan_markers[@]}"; do
        path="$(printf '%s\n' "$marker" | sed -E 's/^<!-- plan-path path=([^ ]+) .*/\1/')"
        class="$(printf '%s\n' "$marker" | sed -E 's/.* class=([^ ]+) destination=.*/\1/')"
        test -z "${plan_classes[$path]+x}"
        plan_classes["$path"]="$class"
    done

    for item in "${actionable[@]}"; do
        path="${item%%$'\t'*}"
        class="${item#*$'\t'}"
        test "${plan_classes[$path]+x}" = x
        test "${plan_classes[$path]}" = "$class"
    done
fi

git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish Residual Review

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The residual inventory and both reports passed post-review verification. Report the report paths, classification totals, unresolved decisions, ordered plan, and developer-confirmed report commit—or `none`—to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The residual review cannot continue safely. Report the exact path, API, ownership, or verification blocker; preserved working-tree state; and maintainer input required.
