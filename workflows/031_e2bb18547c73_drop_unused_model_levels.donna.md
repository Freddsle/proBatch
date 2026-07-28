# Drop Unused Core Model Levels

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Port only the core correction hunks from mixed commit `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`, excluding BERT implementation at hunk level.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-029-e2bb1854-drop-unused-model-levels` and the complete ownership map.
2. Inspect source `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`, parent `20e76c9b9b28c0ec98faa63c3f9382c1347301b9`, only with pinned full-ID Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_effects.R`, `tests/testthat/test-bert-covariates.R`, and `tests/testthat/test-correct_batch_effect.R`.
4. Separate ownership exactly:
   - Companion: BERT helpers drop unused batch levels and expand categorical covariates with at least three levels to `k-1` dummies while keeping binary coding in one column; exclude `R/correct_batch_BERT_PLSDA.R` and `test-bert-covariates.R`.
   - Core: `run_ComBat_core` and `.removeBatchEffect_matrix_step` drop unused batch and covariate factor levels before model construction, transitively affecting `correct_with_ComBat*` and `correct_with_removeBatchEffect*`.
5. Core focused test is `ComBat / removeBatchEffect tolerate unused batch factor levels` in `test-correct_batch_effect.R`. BERT tests remain companion-only.
6. Dependency questions: core uses existing sva, limma, stats, and BiocParallel only; BERT must remain absent from core metadata.
7. Skip only when destination behavior and a focused unused-level regression are already equivalent. Do not skip merely because the source commit is mixed.
8. Never inspect/write external working trees or access `man/**`/`NAMESPACE`. Query Depmesh before every destination edit.
9. If ownership is clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the Core Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare current `R/correct_batch_effects.R::run_ComBat_core` and `.removeBatchEffect_matrix_step` with the exact core hunks.
2. Port only safe unused-level handling and focused core coverage. Preserve current missing-policy, alignment, serial backend, covariate validation, and public signatures.
3. Exclude all BERT helpers, dummy expansion, provider tests, registrations, Roxygen exports, and dependencies.
4. Query Depmesh for `R/correct_batch_effects.R` and `tests/testthat/test-correct_batch_effect.R`; use `apply_patch`.
5. Preserve all user changes. Never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
6. Run the focused correction test when available, parse changed R/test sources, and read `DESCRIPTION`.
7. Create `.session/donna/031-e2bb18547c73-notes.md` with hunk ownership, dependency decisions, checks, and exactly one of each:
   - `<!-- bec-port id=post-sync-029-e2bb1854-drop-unused-model-levels status=changed -->` or `<!-- bec-port id=post-sync-029-e2bb1854-drop-unused-model-levels status=no-change -->`
   - `<!-- focused-check id=post-sync-029-e2bb1854-drop-unused-model-levels source-parse=passed test-parse=passed behavior=passed -->`
8. Then `{{ donna.lib.goto("verify_bec_port") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify the Core Port

```toml donna
id = "verify_bec_port"
kind = "donna.lib.run_script"
save_stdout_to = "bec_verify_stdout"
save_stderr_to = "bec_verify_stderr"
goto_on_success = "review_bec_port"
goto_on_failure = "repair_bec_port"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
notes=".session/donna/031-e2bb18547c73-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port id=post-sync-029-e2bb1854-drop-unused-model-levels status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check id=post-sync-029-e2bb1854-drop-unused-model-levels source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the Core Port

```toml donna
id = "repair_bec_port"
kind = "donna.lib.request_action"
```

Review `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only the core correction hunk, focused core test, or notes after querying Depmesh; preserve BERT exclusions; then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Core Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report core and companion hunk evidence separately, destination changes/no-change, dependency impact, focused verification, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- source-review id=post-sync-029-e2bb1854-drop-unused-model-levels commit=<40-hex-or-none> -->`.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; if corrections are requested, `{{ donna.lib.goto("repair_bec_port") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Reverify the Core Port

```toml donna
id = "reverify_bec_port"
kind = "donna.lib.run_script"
save_stdout_to = "bec_reverify_stdout"
save_stderr_to = "bec_reverify_stderr"
goto_on_success = "compare_split_implementation"
goto_on_failure = "repair_bec_port"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
notes=".session/donna/031-e2bb18547c73-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review id=post-sync-029-e2bb1854-drop-unused-model-levels commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Compare the Split Implementation

```toml donna
id = "compare_split_implementation"
kind = "donna.lib.request_action"
```

1. Use pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split` to read provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`.
2. Core split targets are `R/correct_batch_effects.R::run_ComBat_core`, `.removeBatchEffect_matrix_step`, and `tests/testthat/test-correct_batch_effect.R`.
3. Confirm whether the comparator lacks `droplevels` for batch/covariate factors and lacks the unused-level regression; reassess the source hunk independently.
4. Boundary targets for excluded BERT behavior are `R/registry.R`, `tests/testthat/test-symbol-ownership.R`, and `DESCRIPTION`; confirm BERT implementation/dependency stays absent.
5. Treat the split as rejected evidence and preserve stronger current correction contracts.
6. Never inspect/write external working trees or access generated files.
7. If sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before edits. Apply only independently justified unused-level hardening/focused core tests and preserve existing behavior.
2. Keep all BERT implementation, tests, aliases, and dependencies excluded.
3. Preserve user changes; never stage, commit, amend, or touch generated files.
4. Update notes with per-target rationale and exactly one of each:
   - `<!-- split-stage id=post-sync-029-e2bb1854-drop-unused-model-levels status=changed -->` or `<!-- split-stage id=post-sync-029-e2bb1854-drop-unused-model-levels status=no-change -->`
   - `<!-- split-focused-check id=post-sync-029-e2bb1854-drop-unused-model-levels source-parse=passed test-parse=passed behavior=passed -->`
5. Then `{{ donna.lib.goto("verify_split_adjustment") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify the Split Adjustment

```toml donna
id = "verify_split_adjustment"
kind = "donna.lib.run_script"
save_stdout_to = "split_verify_stdout"
save_stderr_to = "split_verify_stderr"
goto_on_success = "review_split_adjustment"
goto_on_failure = "repair_split_adjustment"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
notes=".session/donna/031-e2bb18547c73-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-stage id=post-sync-029-e2bb1854-drop-unused-model-levels status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check id=post-sync-029-e2bb1854-drop-unused-model-levels source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the Split Adjustment

```toml donna
id = "repair_split_adjustment"
kind = "donna.lib.request_action"
```

Review `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only core correction hardening/tests or notes after querying Depmesh, preserving BERT exclusions; then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split comparison, adopted/rejected core hardening, BERT exclusions, checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- split-review id=post-sync-029-e2bb1854-drop-unused-model-levels commit=<40-hex-or-none> -->`.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; if corrections are requested, `{{ donna.lib.goto("repair_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Reverify the Split Adjustment

```toml donna
id = "reverify_split_adjustment"
kind = "donna.lib.run_script"
save_stdout_to = "split_reverify_stdout"
save_stderr_to = "split_reverify_stderr"
goto_on_success = "finish"
goto_on_failure = "repair_split_adjustment"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
notes=".session/donna/031-e2bb18547c73-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-review id=post-sync-029-e2bb1854-drop-unused-model-levels commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Both review gates passed. Report the core unused-level result, BERT exclusions, and confirmed commit values to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved ownership, dependency, correction, or verification issue and required maintainer input.
