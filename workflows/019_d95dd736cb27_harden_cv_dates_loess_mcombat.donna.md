# Migrate CV, Date, and LOESS Safeguards

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Port the core-owned safeguards from mixed commit `d95dd736cb27d68fb2e21b20ab97d0f69826a663`, exclude mComBat, and compare the pinned split under two maintainer-controlled pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-017-d95dd736-harden-cv-dates-loess-mcombat` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `d95dd736cb27d68fb2e21b20ab97d0f69826a663` against parent `b8a262b4256966d60e1e8452ebde7a1bf471b4af` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/CV_calculation.R`, `R/M_ComBat.R`, `R/date_conversion.R`, and `R/fit_non_linear.R`.
4. Core hunks: `compute_cv` returns `NA` for non-finite/near-zero means; `dates_to_posix` restores locale on exit; `date_to_sample_order` uses first-tie ranking; `loess_regression` and `_opt` suppress extrapolation and return NA fallbacks.
5. Exclude the companion `correct_with_mComBat` sample-count/variance safeguards and all mComBat tests/dependencies.
6. Historical tests changed none. Add focused cases in `test-CV_calculation.R`, `test-date_conversion.R`, and `test-fit_non_linear.R`; no core dependency change is expected.
7. Split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` lacks all four core safeguards in `R/CV_calculation.R`, `R/date_conversion.R`, and `R/fit_non_linear.R`; mComBat is intentionally absent.
8. Query Depmesh for each source/test before editing. Never inspect external working trees, write externally, cherry-pick, or access `man/**` or `NAMESPACE`.
9. If behavior and ownership are clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare current implementations and tests before editing. Port the four independent core safeguards in source order; preserve any newer compatible behavior.
2. Add focused tests for zero/near-zero CV means, locale restoration on success/error, deterministic tied ordering, and LOESS predictions outside fitted support plus warning/error fallback.
3. Confirm no mComBat code or dependency enters core. Use `apply_patch`, preserve user changes, and never cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Create `.session/donna/019-harden-cv-dates-loess-notes.md` with source evidence, Depmesh/dependency results, companion exclusion, checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat -->`
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
notes=".session/donna/019-harden-cv-dates-loess-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat -->' "$notes")" -eq 1
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

Repair only the four core safeguards, focused tests, or notes after querying Depmesh. Preserve the mComBat exclusion. Then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report each safeguard, companion exclusion, tests, dependency result, and working-tree status.
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
notes=".session/donna/019-harden-cv-dates-loess-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat -->' "$notes")" -eq 1
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

1. Compare accepted behavior with commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split`.
2. Inspect the exact core functions and focused tests in `R/CV_calculation.R`, `R/date_conversion.R`, `R/fit_non_linear.R`, `test-CV_calculation.R`, `test-date_conversion.R`, and `test-fit_non_linear.R`.
3. Reassess each absent safeguard independently; never import `R/M_ComBat.R`, inspect/write an external working tree, or access generated files.
4. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if evidence conflicts.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh, then apply only justified additional core hardening/tests. Record each comparator target as adopted, equivalent, or rejected, and keep mComBat absent.
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
notes=".session/donna/019-harden-cv-dates-loess-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat -->' "$notes")" -eq 1
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

Repair only these core safeguards/tests or notes after querying Depmesh; keep mComBat excluded. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split evidence, adopted/rejected hardening, focused tests, and working-tree status.
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
notes=".session/donna/019-harden-cv-dates-loess-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-017-d95dd736-harden-cv-dates-loess-mcombat -->' "$notes")" -eq 1
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

Report the four core safeguard outcomes, mComBat exclusion, and review markers from `.session/donna/019-harden-cv-dates-loess-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved behavior/ownership or verification decision, preserved working tree, and required maintainer input.
