# Migrate removeBatchEffect Missing-Policy Semantics

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Translate core commit `4285c42f31670d2f750dc8eb8c7ff1d0134a342d` into the destination's canonical missing-value policy and compare pinned split hardening under two maintainer-controlled pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `4285c42f31670d2f750dc8eb8c7ff1d0134a342d` against parent `4b117a03f0c5295f2253db26e13b945d7e3027b5` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/correct_batch_effects.R` and `R/handle_missing_values.R`.
4. Core hunks improve removeBatchEffect warnings, make false mean “leave missing values alone,” warn about removals only when they occur, and remove an inaccurate limma claim. Affected functions are `correct_with_removeBatchEffect_dm`, `.removeBatchEffect_matrix_step`, `.run_matrix_method`, and `handle_missing_values`.
5. Historical tests changed none. Focused destination tests are `test-handle_missing_values.R`, `test-correct_batch_effect.R`, `test-batch_effect_steps.R`, and `test-matrix_adapter.R`; there is no dependency change.
6. Compare with split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` targets `R/handle_missing_values.R`, `R/matrix_adapter.R`, and `R/correct_batch_effects.R`. The canonical `error/keep/drop_features/fill` policy plus `fill_value` supersedes legacy boolean syntax.
7. Query every destination artifact through Depmesh before editing. Never inspect external working trees, write externally, cherry-pick, or access `man/**` or `NAMESPACE`.
8. If the semantic mapping is clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare observable warnings and missing-data outcomes, not legacy syntax. Record which intent is already represented by canonical `keep`, `drop_features`, and `fill`.
2. Port only missing core semantics, Roxygen source, and focused tests. Do not reintroduce the boolean interface when the canonical contract covers it.
3. Run narrow missing/correction/adapter tests. Use `apply_patch`, preserve user changes, and never cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Create `.session/donna/016-removebatcheffect-missing-policy-notes.md` with evidence, Depmesh results, semantic/skip decisions, checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy -->`
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
5. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Verify the BEC Port

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
notes=".session/donna/016-removebatcheffect-missing-policy-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy -->' "$notes")" -eq 1
test "$(grep -Ec '^<!-- bec-port status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Fxc '<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the BEC Port

```toml donna
id = "repair_bec_port"
kind = "donna.lib.request_action"
```

Verification failed:

```text
{{ donna.lib.task_variable("bec_verify_stdout") }}
{{ donna.lib.task_variable("bec_verify_stderr") }}
```

Repair only canonical missing-policy behavior, focused tests, Roxygen source, or notes after querying Depmesh. Then `{{ donna.lib.goto("verify_bec_port") }}` or, if an API decision is unresolved, `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report the semantic translation, changes/no-change, focused tests, dependency result, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
3. Do not complete this request until the maintainer explicitly resumes and supplies a 40-hex destination commit or an explicit no-new-commit decision.
4. On resume, add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; for corrections, `{{ donna.lib.goto("repair_bec_port") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Reverify the BEC Port

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
notes=".session/donna/016-removebatcheffect-missing-policy-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy -->' "$notes")" -eq 1
test "$(grep -Ec '^<!-- source-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Compare the Split Implementation

```toml donna
id = "compare_split_implementation"
kind = "donna.lib.request_action"
```

1. Compare the accepted source outcome with split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` through pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split` only.
2. Review canonical policy normalization, `fill_value`, removeBatchEffect's `keep` behavior, removal reporting, and tests across `R/handle_missing_values.R`, `R/matrix_adapter.R`, and `R/correct_batch_effects.R`.
3. Treat the split as rejected evidence; do not copy legacy booleans or inspect/write an external working tree or generated files.
4. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if semantics conflict.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh, then apply only independently justified canonical-policy hardening, source docs, and tests. Record adopted, equivalent, superseded, and rejected behavior.
2. Preserve user changes; use `apply_patch`; never stage, commit, amend, cherry-pick, or touch generated files.
3. Add exactly one of each:
   - `<!-- split-stage status=changed -->` or `<!-- split-stage status=no-change -->`
   - `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`
4. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/016-removebatcheffect-missing-policy-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy -->' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-stage status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Fxc '<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->' "$notes")" -eq 1
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

Split verification failed:

```text
{{ donna.lib.task_variable("split_verify_stdout") }}
{{ donna.lib.task_variable("split_verify_stderr") }}
```

Repair only this workflow's canonical-policy adjustment, tests, source docs, or notes after querying Depmesh. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report the split comparison, adopted/rejected semantics, tests, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
3. Do not complete this request until the maintainer explicitly resumes and supplies a 40-hex destination commit or an explicit no-new-commit decision.
4. On resume, add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; for corrections, `{{ donna.lib.goto("repair_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/016-removebatcheffect-missing-policy-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-014-4285c42f-clarify-removebatcheffect-missing-policy -->' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Report canonical-policy decisions and both maintainer review markers from `.session/donna/016-removebatcheffect-missing-policy-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved missing-policy/API or verification decision, preserved working tree, and required maintainer input.
