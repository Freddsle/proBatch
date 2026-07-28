# Reconcile In-Place Filter Link Effects

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Assess core commit `5edcc4b89f5cbfba07e1ca1057ab475a56501202` without replaying obsolete QFeatures link-removal warnings.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-026-5edcc4b8-report-inplace-link-effects` and the complete ownership map.
2. Inspect source `5edcc4b89f5cbfba07e1ca1057ab475a56501202`, parent `b60e2b169afe45934c18eee1c83469e7e5fed33f`, only with pinned full-ID Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. The sole non-generated path is `R/pb_missing_filters.R`.
4. Core hunks make exported `pb_filterNA` and `pb_groupfilterNA` detect links and warn that in-place filtering removes them; Roxygen describes the effect.
5. No historical focused test was added. Dependencies are existing QFeatures and SummarizedExperiment.
6. Skip rationale: do not replay a warning where current behavior preserves/prunes links. Evaluate the two functions separately and retain only behavior still true.
7. Never inspect/write external working trees or access `man/**`/`NAMESPACE`. Query Depmesh before edits.
8. If clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the Core Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare `pb_filterNA` and `pb_groupfilterNA` separately with current `R/pb_missing_filters.R`; inspect link handling and editable Roxygen sources.
2. Preserve stronger link-pruning behavior. Port a warning or documentation statement only when links are actually removed, and add a focused test for any changed behavior.
3. Query Depmesh before edits; use `apply_patch`.
4. Preserve all user changes. Never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
5. Run narrow missing-helper tests, parse affected files, and read `DESCRIPTION`.
6. Create `.session/donna/028-5edcc4b89f5c-notes.md` with per-function decisions, dependency/test results, and exactly one of each:
   - `<!-- bec-port id=post-sync-026-5edcc4b8-report-inplace-link-effects status=changed -->` or `<!-- bec-port id=post-sync-026-5edcc4b8-report-inplace-link-effects status=no-change -->`
   - `<!-- focused-check id=post-sync-026-5edcc4b8-report-inplace-link-effects source-parse=passed test-parse=passed behavior=passed -->`
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
notes=".session/donna/028-5edcc4b89f5c-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port id=post-sync-026-5edcc4b8-report-inplace-link-effects status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check id=post-sync-026-5edcc4b8-report-inplace-link-effects source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
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

Review `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only link handling/docs, focused tests, or notes after querying Depmesh; then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Core Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report per-function evidence, adopted/skipped behavior, dependency/test impact, checks, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- source-review id=post-sync-026-5edcc4b8-report-inplace-link-effects commit=<40-hex-or-none> -->`.
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
notes=".session/donna/028-5edcc4b89f5c-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review id=post-sync-026-5edcc4b8-report-inplace-link-effects commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
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
2. Split targets are `R/pb_missing_filters.R` and `tests/testthat/test-pb_missing_helpers.R`.
3. Verify that pinned `pb_filterNA` preserves/prunes valid links and tests that behavior; do not replay its obsolete removal warning. Separately assess the remaining `pb_groupfilterNA` removal-warning path and missing focused link coverage.
4. Confirm QFeatures/SummarizedExperiment remain existing core dependencies with no expansion.
5. Never inspect/write external working trees or access generated files.
6. If sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before edits. Adopt only independently justified current link behavior, accurate Roxygen, and focused tests; no-change is valid.
2. Preserve user changes; never stage, commit, amend, or touch generated files.
3. Update notes with separate `pb_filterNA`/`pb_groupfilterNA` rationale and exactly one of each:
   - `<!-- split-stage id=post-sync-026-5edcc4b8-report-inplace-link-effects status=changed -->` or `<!-- split-stage id=post-sync-026-5edcc4b8-report-inplace-link-effects status=no-change -->`
   - `<!-- split-focused-check id=post-sync-026-5edcc4b8-report-inplace-link-effects source-parse=passed test-parse=passed behavior=passed -->`
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
notes=".session/donna/028-5edcc4b89f5c-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-stage id=post-sync-026-5edcc4b8-report-inplace-link-effects status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check id=post-sync-026-5edcc4b8-report-inplace-link-effects source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
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

Review `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only link behavior/docs/tests or notes after querying Depmesh; then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report comparison, adopted/skipped behavior, tests, checks, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. Only then record exactly one `<!-- split-review id=post-sync-026-5edcc4b8-report-inplace-link-effects commit=<40-hex-or-none> -->`.
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
notes=".session/donna/028-5edcc4b89f5c-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-review id=post-sync-026-5edcc4b8-report-inplace-link-effects commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
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

Both review gates passed. Report the per-filter link decisions and confirmed commit values to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved link, dependency, API, or verification issue and required maintainer input.
