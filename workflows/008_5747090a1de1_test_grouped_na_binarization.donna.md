# Migrate Grouped Heatmap Binarization Tests

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate the core test-only behavior from source commit `5747090a1de1cdde780c6db7a84f988aedb9e8af`, with both maintainer pauses even for a no-change result.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-006-5747090a-test-grouped-na-binarization` and the ownership map.
2. Inspect pinned parent `12d4597e9e86ea2daddc77c4d236e692bbfa8fa9`, commit `5747090a1de1cdde780c6db7a84f988aedb9e8af`, and their diff in `/home/yuliya/repos/other/proBatch`.
3. The only changed non-generated path is `tests/testthat/test-plot_missing.R`.
4. Core test hunks cover `.pb_group_missing_matrix` and `plot_grouped_NA_heatmap.default`: `TRUE` means threshold `0.5`, numeric thresholding precedes `drop_complete`, and values outside `[0,1]` fail.
5. Dependency question: none. Ordering question: the historical tests precede implementation commit `7feb6d5cea9ec92d378453ff65d672f231fa44b8`; do not add failing tests if the destination behavior is absent.
6. Exact split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; targets are `R/plot_missing.R` and `tests/testthat/test-plot_missing.R`.
7. Skip rationale: make no R implementation change in this test-only source stage; record no-change if equivalent tests already exist or if implementation must be handled by the next workflow.
8. Query Depmesh before edits. Use only full-ID Git object reads externally; never inspect/write external working trees or access `man/**`/`NAMESPACE`. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
9. Continue to `{{ donna.lib.goto("port_bec_behavior") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

Port only safe focused tests after confirming the implementation contract. Preserve user changes; never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files. Create `.session/donna/008-5747090a-test-grouped-na-binarization-notes.md` with evidence, ordering/dependency decision, Depmesh queries, checks, and exactly one marker of each form:

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
notes=".session/donna/008-5747090a-test-grouped-na-binarization-notes.md"
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

Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only the focused tests or notes after Depmesh queries. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. Report tests, ordering decision, checks, and status. The agent must not stage, commit, amend, or draft a commit message. Do not complete until the maintainer explicitly resumes and confirms a commit SHA or no-new-commit decision. Add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_bec_port") }}`, `{{ donna.lib.goto("repair_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/008-5747090a-test-grouped-na-binarization-notes.md"
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

Use `git -C /home/yuliya/repos/other/proBatch-core-split show`, `diff`, `ls-tree`, or `grep` only against provenance object `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`. Confirm the implementation exists and the exact tests are absent. Never inspect/write the external working tree or generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

After Depmesh queries, add only independently justified tests, or record no-change. Preserve user changes; never stage, commit, amend, or touch generated files. Add exactly one concrete split-stage status marker and `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/008-5747090a-test-grouped-na-binarization-notes.md"
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

Inspect `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only these tests/notes after Depmesh queries. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause, including no-change. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_split_adjustment") }}`, `{{ donna.lib.goto("repair_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/008-5747090a-test-grouped-na-binarization-notes.md"
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

Report the test-only decision, ordering rationale, and both review markers.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved test/implementation ordering or verification decision.
