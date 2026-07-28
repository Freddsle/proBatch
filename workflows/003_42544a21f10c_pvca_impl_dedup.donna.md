# Migrate PVCA Implementation Consolidation

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate the core-owned behavior from source commit `42544a21f10ca6960d3e4c44d2833f764054d721`, then independently compare the stopped split under two maintainer-controlled review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-001-42544a21-pvca-impl-dedup` and the complete ownership map.
2. Inspect only pinned Git objects in `/home/yuliya/repos/other/proBatch`: parent `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`, commit `42544a21f10ca6960d3e4c44d2833f764054d721`, and their diff.
3. Changed non-generated paths are `R/explained_variance_plots.R` and `R/proteome_wide_diagnostics.R`.
4. Core hunks capture six PVCA method bodies as `.pb_*_impl` aliases, delegate duplicate definitions, and make `plot_PVCA.df` an S3 generic. Inspect `calculate_PVCA.default`, `calculate_PVCA.ProBatchFeatures`, `plot_PVCA.default`, `plot_PVCA.ProBatchFeatures`, `prepare_PVCA_df.default`, `prepare_PVCA_df.ProBatchFeatures`, and `plot_PVCA.df`.
5. No historical test changed. Compare focused PVCA coverage in `tests/testthat/test-proteome_wide_diagnostics.R` and `tests/testthat/test-symbol-ownership.R`. Confirm whether existing `pvca`, `Biobase`, `ggplot2`, `dplyr`, and `gridExtra` declarations suffice.
6. Split comparator is exactly `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; targets are `R/proteome_wide_diagnostics.R` and `tests/testthat/test-symbol-ownership.R`.
7. Skip rationale: do not recreate `R/explained_variance_plots.R` or the alias mechanism when the destination already has one authoritative implementation; migrate semantic behavior and coverage only.
8. Query all Depmesh relations for every destination artifact before editing.
9. External reads must use `git -C ... show`, `diff`, `ls-tree`, or `grep` against full commit IDs. Never inspect or write an external working tree. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. Never read or enumerate `man/**` or `NAMESPACE`. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, or touch generated files. Draft a commit message only at a maintainer review gate as required below.
10. If ownership and behavior are clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare the retained PVCA semantics with the current destination before editing.
2. Port only missing core behavior, Roxygen sources under `R/`, and focused tests. Preserve all user changes.
3. Do not cherry-pick, stage, commit, amend, draft a commit message, or touch generated files. Use `apply_patch`.
4. Run narrow agent-safe checks and create `.session/donna/003-42544a21-pvca-impl-dedup-notes.md` with evidence, Depmesh queries, decisions, commands/results, and exactly one of each marker:
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
5. Record the consolidation-based skip rationale for every no-change decision.
6. When ready, `{{ donna.lib.goto("verify_bec_port") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/003-42544a21-pvca-impl-dedup-notes.md"
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

Verification failed. Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only this port and its notes after any required Depmesh query. Preserve ownership exclusions and user changes. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}` if judgment is unresolved.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is a hard maintainer review and commit pause.

1. Report source evidence, changed/no-change outcome, dependencies, focused checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action request until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. On resume, add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker to the notes.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; if corrections are requested, `{{ donna.lib.goto("repair_bec_port") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/003-42544a21-pvca-impl-dedup-notes.md"
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

1. Read the provenance note only as pinned object `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` in `/home/yuliya/repos/other/proBatch-core-split`.
2. Compare with exact split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db`, focusing on `R/proteome_wide_diagnostics.R`, PVCA tests, and `tests/testthat/test-symbol-ownership.R`.
3. Treat consolidation as evidence, not an accepted patch. Never inspect generated files or write externally.
4. If evidence is sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh, then apply only independently justified PVCA hardening or focused tests.
2. Preserve user changes; do not stage, commit, amend, or touch generated files.
3. Update the notes with rationale and exactly one of each:
   - `<!-- split-stage status=changed -->` or `<!-- split-stage status=no-change -->`
   - `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`
4. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/003-42544a21-pvca-impl-dedup-notes.md"
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

Inspect `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only justified split-stage work after Depmesh queries. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is a hard maintainer review and commit pause.

1. Report comparator evidence, adopted/rejected adjustments, checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action request until the maintainer explicitly resumes and confirms a destination commit SHA or explicit no-new-commit decision.
4. On resume, add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; if corrections are requested, `{{ donna.lib.goto("repair_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/003-42544a21-pvca-impl-dedup-notes.md"
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

Both maintainer reviews and reverification gates passed. Report source and split decisions and the confirmed commit markers to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved PVCA ownership, API, dependency, or verification decision, preserved user changes, and exact maintainer input required.
