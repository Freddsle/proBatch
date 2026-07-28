# Migrate Grouped Density Color Schemes

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate grouped-density color behavior from source commit `0dbce8129444c9079a560fb9bcd60b328b1f054c` and compare the stopped split under two maintainer pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-008-0dbce812-density-color-scheme` and the ownership map.
2. Inspect pinned parent `7feb6d5cea9ec92d378453ff65d672f231fa44b8`, commit `0dbce8129444c9079a560fb9bcd60b328b1f054c`, and diff in `/home/yuliya/repos/other/proBatch`.
3. Changed non-generated paths: `R/plot_missing.R` and `tests/testthat/test-plot_missing.R`.
4. Core hunks add brewer/named/list `color_scheme`, raw group keys, a common theme, and complete Roxygen parameters for `plot_NA_density.default`, `plot_NA_density.ProBatchFeatures`, `.pb_plot_missing_density`, `.pb_prepare_grouped_missing_annotation`, `.pb_resolve_missing_density_color_scheme`, and `.pb_missing_density_theme`.
5. Focused tests update group-key expectations and verify named colors with `ggplot_build`.
6. Dependency question: use existing `ggplot2` and core color helpers; no package metadata change should be needed.
7. Exact split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; target `R/plot_missing.R` contains internal helpers but lacks the public grouped-density signature and test.
8. Skip rationale: helper presence alone is not sufficient; skip duplicate/dead helper edits and port only coherent public behavior and coverage if accepted.
9. Query Depmesh before edits. External reads are pinned Git objects only; never inspect/write external working trees or `man/**`/`NAMESPACE`. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, or touch generated files. Draft a commit message only at a maintainer review gate as required below.
10. Continue to `{{ donna.lib.goto("port_bec_behavior") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

Compare and port only coherent missing core API/Roxygen/tests. Preserve user changes; no cherry-pick, staging, commit, amend, commit-message drafting, or generated files. Create `.session/donna/010-0dbce812-density-color-scheme-notes.md` with evidence, dependency/API decision, Depmesh queries, checks, and exactly one:

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
notes=".session/donna/010-0dbce812-density-color-scheme-notes.md"
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

Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only this API/test/notes scope after Depmesh queries. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_bec_port") }}`, `{{ donna.lib.goto("repair_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/010-0dbce812-density-color-scheme-notes.md"
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

Use `git -C /home/yuliya/repos/other/proBatch-core-split show`, `diff`, `ls-tree`, or `grep` only against provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and exact comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`. Compare public signatures, helper reachability, Roxygen sources, and focused tests. Never write externally or inspect generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

After Depmesh queries, apply only coherent and independently justified color-scheme behavior/tests. Preserve user changes; no staging, commit, amend, or generated files. Add one concrete split-stage marker and `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/010-0dbce812-density-color-scheme-notes.md"
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

Hard maintainer review and commit pause. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_split_adjustment") }}`, `{{ donna.lib.goto("repair_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/010-0dbce812-density-color-scheme-notes.md"
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

Report public color-scheme and helper-reachability decisions and both review markers.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved API, dependency, or verification decision.
