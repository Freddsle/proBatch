# Migrate All-NA Helper Documentation Metadata

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Audit documentation-only commit `77a3365f0e13e9abca40a14e1898a351a1cdbc7f` and preserve its intent only if the preceding internal-helper design survives, with two maintainer-controlled review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-012-77a3365f-demote-all-na-helper-roxygen` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `77a3365f0e13e9abca40a14e1898a351a1cdbc7f` against parent `fdafc6c8fe1391de1cc0d2db36d346cd0d4e61d5` using only pinned `git -C /home/yuliya/repos/other/proBatch show`, `diff`, `ls-tree`, or `grep`.
3. The sole non-generated path is `R/utility_funcs.R`. The hunk converts Roxygen blocks for `.pb_strip_allna` and `.pb_restore_allna` to ordinary comments; it changes no runtime behavior, export, test, or dependency.
4. The pinned split comparator is `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; it has neither helper, so there is no current split implementation target. Apply this metadata only if workflow 013 adopted those helpers.
5. Query all Depmesh relations for `R/utility_funcs.R` before editing. Never inspect external working trees, write externally, cherry-pick, or access `man/**` or `NAMESPACE`.
6. If the prerequisite and skip decision are clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Check whether the two helpers exist in the accepted destination source. If absent, record `no-change`; if present, preserve their internal ordinary-comment status without generating documentation.
2. No focused runtime test or package dependency change is expected. Parse any edited R source and perform a source-level documentation check only.
3. Use `apply_patch`, preserve user changes, and do not cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Create `.session/donna/014-demote-all-na-helper-roxygen-notes.md` with prerequisite evidence, Depmesh output, the skip/apply decision, checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-012-77a3365f-demote-all-na-helper-roxygen -->`
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
5. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}` if unresolved.

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
notes=".session/donna/014-demote-all-na-helper-roxygen-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-012-77a3365f-demote-all-na-helper-roxygen -->' "$notes")" -eq 1
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

Repair only this metadata decision or its notes after any required Depmesh query. Preserve the no-generated-files rule. Then `{{ donna.lib.goto("verify_bec_port") }}` or, if unresolved, `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report the prerequisite result, changed/no-change decision, source check, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
3. Do not complete this request until the maintainer explicitly resumes and supplies a 40-hex destination commit or an explicit no-new-commit decision.
4. On resume, add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker.
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
notes=".session/donna/014-demote-all-na-helper-roxygen-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-012-77a3365f-demote-all-na-helper-roxygen -->' "$notes")" -eq 1
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

Use only pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split` to confirm that split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` has no helper or documentation target. Record the resulting no-change rationale unless the accepted source stage introduced the helpers. Never inspect/write the split working tree or generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if evidence conflicts.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before any edit. Normally record no split change; if helpers now exist, ensure they remain internal and source-documented without generated output.
2. Preserve user changes and never stage, commit, amend, cherry-pick, or touch generated files.
3. Add exactly one of each to the notes:
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
notes=".session/donna/014-demote-all-na-helper-roxygen-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-012-77a3365f-demote-all-na-helper-roxygen -->' "$notes")" -eq 1
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

Repair only this workflow's metadata/notes after any required Depmesh query. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report the split no-target evidence, any adjustment, checks, and working-tree status.
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
notes=".session/donna/014-demote-all-na-helper-roxygen-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-012-77a3365f-demote-all-na-helper-roxygen -->' "$notes")" -eq 1
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

Report the documentation metadata decision and both maintainer-confirmed review markers from `.session/donna/014-demote-all-na-helper-roxygen-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved prerequisite or verification decision, preserved working tree, and exact maintainer input required.
