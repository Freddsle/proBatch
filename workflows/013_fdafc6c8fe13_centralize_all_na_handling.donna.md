# Migrate Centralized All-NA Handling

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Audit mixed commit `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5`, port only justified core behavior, and compare it with split-core commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` under two maintainer-controlled review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-011-fdafc6c8-centralize-all-na-handling` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5` against parent `7db05faed18f75572d92553ae668288c326bda2a` in `/home/yuliya/repos/other/proBatch`. External repository reads are limited to `git -C ... show`, `diff`, `ls-tree`, or `grep` against those full IDs.
3. The non-generated source paths are `R/ProBatchFeatures.R`, `R/correct_batch_RUVIIIC.R`, `R/correct_batch_effects.R`, `R/fit_non_linear.R`, `R/impute_missForest.R`, `R/impute_omicsGMF.R`, `R/normalize.R`, and `R/utility_funcs.R`.
4. Retain for review only the core hunks for `.pb_strip_allna`, `.pb_restore_allna`, registry median-normalization no-handling semantics, `center_feature_batch` wide validation, `correct_batch_effects_df/dm` forwarding, nullable group handling, and the `fit_nonlinear` comment. Exclude at hunk level the RUVIIIC, missForest, and omicsGMF consumers and their provider dependencies.
5. Focused targets are `tests/testthat/test-correct_batch_effect.R`, `test-batch_effect_steps.R`, `test-handle_missing_values.R`, `test-matrix_adapter.R`, `test-normalize.R`, `test-fit_non_linear.R`, and `test-registry.R`. Split targets are `R/correct_batch_effects.R`, `R/handle_missing_values.R`, `R/matrix_adapter.R`, `R/normalize.R`, and `R/registry.R`.
6. Resolve whether the canonical missing policy supersedes each boolean hunk, whether helper code has any surviving core consumer, and whether wrapper forwarding remains compatible. Provider dependencies must remain absent from core.
7. Never inspect an external working tree, write externally, cherry-pick, or access `man/**` or `NAMESPACE`. Query all Depmesh relations for each destination artifact before editing.
8. If ownership and compatibility are resolved, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare each retained hunk with the current destination and record already-equivalent or superseded behavior before editing. Do not recreate provider-only helpers without a core consumer.
2. Port only missing core behavior, editable Roxygen source, and focused tests. Use `apply_patch`, preserve all user changes, and never copy a mixed source file wholesale.
3. Do not cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Run narrow agent-safe checks and create `.session/donna/013-centralize-all-na-handling-notes.md` containing evidence, Depmesh results, dependency and skip decisions, commands/results, and exactly one of each:
   - `<!-- migration-id id=post-sync-011-fdafc6c8-centralize-all-na-handling -->`
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
5. When ready, `{{ donna.lib.goto("verify_bec_port") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/013-centralize-all-na-handling-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-011-fdafc6c8-centralize-all-na-handling -->' "$notes")" -eq 1
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

Verification failed.

```text
{{ donna.lib.task_variable("bec_verify_stdout") }}
{{ donna.lib.task_variable("bec_verify_stderr") }}
```

Repair only this commit's core port, tests, or notes; query Depmesh before additional governed edits and preserve all exclusions. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}` if an ownership/API decision is unresolved.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report source evidence, core/companion decisions, changes or no-change, dependencies, focused checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action request until the maintainer explicitly resumes and supplies a 40-hex destination commit or explicitly chooses no new commit.
4. On explicit resume, add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker to the notes.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; for corrections, `{{ donna.lib.goto("repair_bec_port") }}`; if rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/013-centralize-all-na-handling-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-011-fdafc6c8-centralize-all-na-handling -->' "$notes")" -eq 1
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

1. Compare the accepted source-stage outcome with immutable split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` using only pinned `git -C /home/yuliya/repos/other/proBatch-core-split show`, `diff`, `ls-tree`, or `grep`.
2. Confirm that wide validation is already present, canonical missing handling supersedes legacy booleans, all-NA helpers have no retained consumer, group-aware normalization was removed, and wrapper forwarding still requires an explicit API decision.
3. Treat the split as rejected evidence; never inspect/write its working tree or access generated files.
4. If evidence is sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh for every destination artifact, then apply only independently justified core hardening and tests. Record adopted, equivalent, superseded, deferred, and companion-excluded targets.
2. Preserve user changes; use `apply_patch`; never stage, commit, amend, cherry-pick, or touch generated files.
3. Update the notes with exactly one of each:
   - `<!-- split-stage status=changed -->` or `<!-- split-stage status=no-change -->`
   - `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`
4. When ready, `{{ donna.lib.goto("verify_split_adjustment") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/013-centralize-all-na-handling-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-011-fdafc6c8-centralize-all-na-handling -->' "$notes")" -eq 1
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

Split verification failed.

```text
{{ donna.lib.task_variable("split_verify_stdout") }}
{{ donna.lib.task_variable("split_verify_stderr") }}
```

Repair only this workflow's justified split changes, tests, or notes; query Depmesh before more edits and retain all exclusions. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split evidence, adopted/rejected hardening, changes or no-change, checks, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
3. Do not complete this action request until the maintainer explicitly resumes and supplies a 40-hex destination commit or explicitly chooses no new commit.
4. On explicit resume, add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker to the notes.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; for corrections, `{{ donna.lib.goto("repair_split_adjustment") }}`; if rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/013-centralize-all-na-handling-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-011-fdafc6c8-centralize-all-na-handling -->' "$notes")" -eq 1
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

Both review gates and reverification passed. Report source/split decisions and maintainer-confirmed commit markers from `.session/donna/013-centralize-all-na-handling-notes.md` to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the exact unresolved API, ownership, dependency, or verification decision, confirm the working tree was preserved, and identify the maintainer input required.
