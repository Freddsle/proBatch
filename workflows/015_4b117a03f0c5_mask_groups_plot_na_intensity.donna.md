# Migrate Group Masking and NA-Intensity Plotting

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate core commit `4b117a03f0c5295f2253db26e13b945d7e3027b5`, recover any missing grouped-filter tests, and compare the implementation with pinned split-core evidence under two maintainer review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-013-4b117a03-mask-groups-plot-na-intensity` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `4b117a03f0c5295f2253db26e13b945d7e3027b5` against parent `77a3365f0e13e9abca40a14e1898a351a1cdbc7f` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/pb_missing_filters.R`, added `R/plot_NA_intensity.R`, `tests/testthat/test-pb_missing_helpers.R`, and added `tests/testthat/test-plot_NA_intensity.R`.
4. Core behavior adds/logs `mask_failing` in `pb_groupfilterNA` and adds `plot_NA_intensity`, its default/`ProBatchFeatures` methods, grouped statistics, spline trend, Spearman labels, and faceting.
5. Existing dependencies are SummarizedExperiment, ggplot2, scales, and splines; determine whether current `DESCRIPTION` already declares everything retained. Focused tests are the two changed historical test files.
6. Split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` already contains both implementations in `R/pb_missing_filters.R` and `R/plot_NA_intensity.R`, plus plot tests, but lacks the two focused `mask_failing` cases.
7. Query all Depmesh relations for destination sources/tests before editing. Never inspect external working trees, write externally, cherry-pick, or access `man/**` or `NAMESPACE`.
8. If evidence is sufficient, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare behavior and signatures with the current destination before editing. Preserve equivalent implementations and port only missing behavior or tests.
2. Ensure focused coverage verifies default masking, `mask_failing = FALSE`, and operation-log recording; run `test-pb_missing_helpers.R` and `test-plot_NA_intensity.R` narrowly when available.
3. Use `apply_patch`, preserve all user changes, and do not cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Create `.session/donna/015-mask-groups-plot-na-intensity-notes.md` with evidence, Depmesh/dependency decisions, checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-013-4b117a03-mask-groups-plot-na-intensity -->`
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
notes=".session/donna/015-mask-groups-plot-na-intensity-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-013-4b117a03-mask-groups-plot-na-intensity -->' "$notes")" -eq 1
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

Repair only group masking, NA-intensity behavior/tests, dependency metadata, or notes after querying Depmesh. Then `{{ donna.lib.goto("verify_bec_port") }}` or, for an unresolved API/dependency issue, `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report source evidence, changes or equivalence, dependency impact, focused tests, and working-tree status.
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
notes=".session/donna/015-mask-groups-plot-na-intensity-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-013-4b117a03-mask-groups-plot-na-intensity -->' "$notes")" -eq 1
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

1. Compare accepted behavior with split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split`.
2. Review `R/pb_missing_filters.R`, `R/plot_NA_intensity.R`, `tests/testthat/test-pb_missing_helpers.R`, and `tests/testthat/test-plot_NA_intensity.R`. Confirm implementation equivalence and specifically locate the absent default/false masking cases.
3. Treat split code as rejected evidence; do not copy wholesale, inspect/write an external working tree, or access generated files.
4. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if evidence conflicts.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh, then add only independently justified missing tests or hardening. Record every split target as adopted, already equivalent, or rejected.
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
notes=".session/donna/015-mask-groups-plot-na-intensity-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-013-4b117a03-mask-groups-plot-na-intensity -->' "$notes")" -eq 1
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

Repair only this workflow's split-stage behavior/tests or notes after any required Depmesh query. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split comparison, adopted/rejected changes, focused tests, and working-tree status.
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
notes=".session/donna/015-mask-groups-plot-na-intensity-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-013-4b117a03-mask-groups-plot-na-intensity -->' "$notes")" -eq 1
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

Report implementation/test decisions and both maintainer review markers from `.session/donna/015-mask-groups-plot-na-intensity-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved behavior, dependency, or verification decision, preserved working tree, and required maintainer input.
