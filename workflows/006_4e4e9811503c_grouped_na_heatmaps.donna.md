# Migrate Grouped Missingness Heatmaps

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate grouped missingness heatmap behavior from source commit `4e4e9811503cc7da1e35a8113a2fb8383fc5007b` and compare the stopped split under two maintainer pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-004-4e4e9811-grouped-na-heatmaps` and the ownership map.
2. Inspect pinned parent `65f70a46c4cf44e2717744aeafbf8acbe83b0378`, commit `4e4e9811503cc7da1e35a8113a2fb8383fc5007b`, and diff in `/home/yuliya/repos/other/proBatch`.
3. Changed non-generated paths: `R/plot_missing.R` and `tests/testthat/test-plot_missing.R`.
4. Core hunks extract heatmap drawing, annotation alignment, and arrangement helpers; add `plot_grouped_NA_heatmap`, default and `ProBatchFeatures` methods; group one or more metadata columns into observed fractions and group-size labels; and refactor `plot_NA_heatmap.default`/`.ProBatchFeatures`.
5. Focused tests: one-/multi-column aggregation, required `color_by`, default grouped plot, and multi-assay arrangement.
6. Dependency question: retain existing `pheatmap`, `gridExtra`, `grid`, `grDevices`, and `SummarizedExperiment` only; determine whether any metadata change is actually needed.
7. Exact split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; targets: `R/plot_missing.R` and `tests/testthat/test-plot_missing.R`.
8. Skip rationale: do not duplicate already-retained grouped API/helpers; port only missing semantics and exact focused coverage.
9. Query Depmesh before edits. External access is pinned Git-object reads only, no external writes or working-tree inspection, and no `man/**`/`NAMESPACE`. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
10. Continue to `{{ donna.lib.goto("port_bec_behavior") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

Compare first; port missing grouped heatmap semantics, Roxygen source, and focused tests only. Preserve user changes. No cherry-pick, stage, commit, amend, commit-message drafting, or generated files; use `apply_patch`.

Create `.session/donna/006-4e4e9811-grouped-na-heatmaps-notes.md` with evidence, Depmesh queries, dependency/no-change rationale, checks, and exactly one:

- `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
- `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`

Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Verify BEC Port

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
notes=".session/donna/006-4e4e9811-grouped-na-heatmaps-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair BEC Port

```toml donna
id = "repair_bec_port"
kind = "donna.lib.request_action"
```

Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only grouped heatmap behavior/tests/notes after Depmesh queries. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. Report behavior, tests, dependencies, status, and no-change rationale where applicable. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Record `<!-- source-review commit=<40-hex-or-none> -->` exactly once. Then `{{ donna.lib.goto("reverify_bec_port") }}`, `{{ donna.lib.goto("repair_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Reverify BEC Port

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
notes=".session/donna/006-4e4e9811-grouped-na-heatmaps-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- source-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Compare Split Implementation

```toml donna
id = "compare_split_implementation"
kind = "donna.lib.request_action"
```

Use only pinned split provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` in `/home/yuliya/repos/other/proBatch-core-split`. Compare grouped API/helpers and the absent focused tests in `R/plot_missing.R` and `tests/testthat/test-plot_missing.R`. Treat split code as rejected evidence; never write externally or inspect generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

Query Depmesh and apply only independently justified grouped-heatmap hardening/tests. Preserve user changes; no staging, commit, amend, or generated files. Add one concrete split-stage status marker and exactly one `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->` marker to the notes. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Verify Split Adjustment

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
notes=".session/donna/006-4e4e9811-grouped-na-heatmaps-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-stage status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair Split Adjustment

```toml donna
id = "repair_split_adjustment"
kind = "donna.lib.request_action"
```

Inspect `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only this adjustment after Depmesh queries. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Record `<!-- split-review commit=<40-hex-or-none> -->` exactly once. Then `{{ donna.lib.goto("reverify_split_adjustment") }}`, `{{ donna.lib.goto("repair_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Reverify Split Adjustment

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
notes=".session/donna/006-4e4e9811-grouped-na-heatmaps-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-stage status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Report grouped-heatmap source/split decisions and both review markers.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved API, dependency, or verification decision and required maintainer input.
