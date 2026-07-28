# Preserve Explicit Final Log Assays

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Assess core commit `b60e2b169afe45934c18eee1c83469e7e5fed33f` and retain its final-log materialization behavior without regressing split lineage hardening.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-025-b60e2b16-store-explicit-final-log-assay` and the complete ownership map.
2. Inspect source `b60e2b169afe45934c18eee1c83469e7e5fed33f`, parent `96d38eb7449e4d38c0d6a3fffd66e17f145f669f`, only with pinned full-ID Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/ProBatchFeatures.R` and `tests/testthat/test-ProBatchFeatures.R`.
4. Core hunks: log/log2 registry functions accept `log_base`, `offset`, and legacy `base`/`pseudo`; exported `pb_transform` stores an otherwise ephemeral final log step when `final_name` is given and does not rename the source assay.
5. Focused test: `pb_transform with log2 and final_name stores transformed (not raw) values`.
6. Dependency question: behavior uses base R and existing QFeatures/SummarizedExperiment/S4Vectors infrastructure; no new dependency is justified.
7. Compare for existing equivalence before editing. Skip/no-change is required when current registry and `pb_transform` already preserve values, operation log, lineage, and final naming.
8. Never inspect/write external working trees or access `man/**`/`NAMESPACE`. Query Depmesh for every destination artifact.
9. If clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the Core Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare source behavior with `R/ProBatchFeatures.R`, current registry implementation, and focused tests. Detect already-equivalent or stronger behavior first.
2. Port only a missing core behavior and its focused test/Roxygen source. Do not overwrite stricter assay-name validation or lineage/log synchronization.
3. Query Depmesh before edits; use `apply_patch`.
4. Preserve all user changes. Never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
5. Run the narrow focused test if available, parse affected R/test files, and read `DESCRIPTION`.
6. Create `.session/donna/027-b60e2b169afe-notes.md` with equivalence evidence, dependency/test results, and exactly one of each:
   - `<!-- bec-port id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=changed -->` or `<!-- bec-port id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=no-change -->`
   - `<!-- focused-check id=post-sync-025-b60e2b16-store-explicit-final-log-assay source-parse=passed test-parse=passed behavior=passed -->`
7. Then `{{ donna.lib.goto("verify_bec_port") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/027-b60e2b169afe-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check id=post-sync-025-b60e2b16-store-explicit-final-log-assay source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
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

Review `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only this core port, focused test, or notes after querying Depmesh; then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Core Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report evidence, change/equivalence decision, dependency impact, focused checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- source-review id=post-sync-025-b60e2b16-store-explicit-final-log-assay commit=<40-hex-or-none> -->`.
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
notes=".session/donna/027-b60e2b169afe-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review id=post-sync-025-b60e2b16-store-explicit-final-log-assay commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
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
2. Split targets are `R/registry.R::.pb_log2_step`, `.pb_log_step`, `.pb_register_core_steps`; `R/ProBatchFeatures.R::pb_transform`; and `tests/testthat/test-ProBatchFeatures.R`.
3. Verify explicit final log assays are materialized with correct transformed values, validated names, operation-log destination, lineage chain, and pipeline name.
4. Treat the split as rejected evidence; do not regress stronger current behavior.
5. Never inspect/write external working trees or access generated files.
6. If sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before edits. Apply only independently justified hardening/tests not already equivalent; no change is expected if the stronger split behavior is present.
2. Preserve user changes; never stage, commit, amend, or touch generated files.
3. Update notes with target-by-target rationale and exactly one of each:
   - `<!-- split-stage id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=changed -->` or `<!-- split-stage id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=no-change -->`
   - `<!-- split-focused-check id=post-sync-025-b60e2b16-store-explicit-final-log-assay source-parse=passed test-parse=passed behavior=passed -->`
4. Then `{{ donna.lib.goto("verify_split_adjustment") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/027-b60e2b169afe-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-stage id=post-sync-025-b60e2b16-store-explicit-final-log-assay status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check id=post-sync-025-b60e2b16-store-explicit-final-log-assay source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
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

Review `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only this adjustment, focused test, or notes after querying Depmesh; then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report comparison, adopted/rejected hardening, checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- split-review id=post-sync-025-b60e2b16-store-explicit-final-log-assay commit=<40-hex-or-none> -->`.
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
notes=".session/donna/027-b60e2b169afe-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-review id=post-sync-025-b60e2b16-store-explicit-final-log-assay commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
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

Both review gates passed. Report the retained/equivalent behavior and confirmed commit values to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved API, dependency, lineage, or verification issue and required maintainer input.
