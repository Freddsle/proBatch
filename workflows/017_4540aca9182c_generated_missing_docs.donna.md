# Record Generated Missing-Documentation No-Change

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Record the generated-only disposition of commit `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57` without inspecting generated artifacts, while preserving both maintainer review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-015-4540aca9-generated-missing-docs` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit metadata for `4540aca9182c6708fe9bda0b8fc33d2cf8c13e57` and parent `4285c42f31670d2f750dc8eb8c7ff1d0134a342d` with pinned `git -C /home/yuliya/repos/other/proBatch show --no-patch`. Confirm the generated-file-excluding diff is empty without enumerating excluded paths.
3. There are no non-generated changed paths, core functions, focused tests, or dependencies. The editable-source intent was already accounted for by workflow 016.
4. Split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` has no editable target. Documentation regeneration remains maintainer-owned.
5. Never read, enumerate, compare, validate, or migrate `man/**` or `NAMESPACE`; never inspect/write an external working tree, cherry-pick, stage, commit, or amend.
6. If the generated-only evidence is confirmed, `{{ donna.lib.goto("port_bec_behavior") }}`; if it conflicts, `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Make no package-source change. Record why generated output is skipped and verify that workflow 016 owns the editable-source decision.
2. Create `.session/donna/017-generated-missing-docs-notes.md` with evidence and exactly:
   - `<!-- migration-id id=post-sync-015-4540aca9-generated-missing-docs -->`
   - `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
3. Preserve all user changes. Do not stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}` if the no-change disposition cannot be justified.

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
notes=".session/donna/017-generated-missing-docs-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-015-4540aca9-generated-missing-docs -->' "$notes")" -eq 1
test "$(grep -Fxc '<!-- bec-port status=no-change -->' "$notes")" -eq 1
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

Repair only the no-change evidence or session markers. Do not inspect or modify generated artifacts. Then `{{ donna.lib.goto("verify_bec_port") }}` or, if evidence conflicts, `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report the generated-only evidence, no-change result, verification, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
3. Do not complete this request until the maintainer explicitly resumes and confirms an explicit no-new-commit decision or supplies a 40-hex destination commit.
4. On resume, add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; for note corrections, `{{ donna.lib.goto("repair_bec_port") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/017-generated-missing-docs-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-015-4540aca9-generated-missing-docs -->' "$notes")" -eq 1
test "$(grep -Fxc '<!-- bec-port status=no-change -->' "$notes")" -eq 1
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

Confirm via pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split` that commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` provides no editable source target for this generated-only commit. Do not inspect any generated path or external working tree. Record no-change, then `{{ donna.lib.goto("apply_split_adjustment") }}`; if evidence conflicts, `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Make no package-source change and add exactly:
   - `<!-- split-stage status=no-change -->`
   - `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`
2. Preserve user changes and never stage, commit, amend, cherry-pick, regenerate, or touch generated files.
3. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/017-generated-missing-docs-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-015-4540aca9-generated-missing-docs -->' "$notes")" -eq 1
test "$(grep -Fxc '<!-- split-stage status=no-change -->' "$notes")" -eq 1
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

Repair only no-change evidence or markers, without inspecting generated files. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report no-target evidence, no-change verification, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
3. Do not complete this request until the maintainer explicitly resumes and confirms no new commit or supplies a 40-hex destination commit.
4. On resume, add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; for note corrections, `{{ donna.lib.goto("repair_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/017-generated-missing-docs-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-015-4540aca9-generated-missing-docs -->' "$notes")" -eq 1
test "$(grep -Fxc '<!-- split-stage status=no-change -->' "$notes")" -eq 1
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

Report the confirmed generated-only no-change outcome and review markers from `.session/donna/017-generated-missing-docs-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report conflicting evidence or marker failure, confirm generated files and the working tree were preserved, and identify required maintainer input.
