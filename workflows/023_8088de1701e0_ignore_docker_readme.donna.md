# Review Docker Ignore Metadata

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Review metadata commit `8088de1701e0908cca25b978105f8d6c7bfccc20` against the current core and stopped split, with two maintainer-controlled review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-021-8088de17-ignore-docker-readme` and the complete API ownership map.
2. Inspect source commit `8088de1701e0908cca25b978105f8d6c7bfccc20`, parent `72d11d1f7f92a10cb9a5244e68f9c19e830197b9`, using only pinned `git -C /home/yuliya/repos/other/proBatch show`, `diff`, `ls-tree`, or `grep` commands with full object IDs.
3. The complete non-generated path set is `.gitignore`. The hunk merely reorders existing `Dockerfile`/`.dockerignore` ignores and adds `Docker_README.md`; there are no core functions, tests, or package dependencies.
4. Classification is metadata. The expected source-stage result is no change unless current repository policy independently requires `Docker_README.md`; reordering alone is not a port.
5. Never inspect either external working tree, write externally, or read/enumerate `man/**` or `NAMESPACE`.
6. Query all Depmesh relations before editing `.gitignore`.
7. If the evidence is complete, `{{ donna.lib.goto("port_bec_behavior") }}`; if policy cannot be resolved, `{{ donna.lib.goto("blocked") }}`.

## Record the Source Decision

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare the metadata intent with current `.gitignore`; prefer an explicit no-change decision when the extra ignore is obsolete or unnecessary.
2. If independently justified, edit only `.gitignore` with `apply_patch`.
3. Preserve all user changes. Never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
4. Read `DESCRIPTION` to confirm there is no dependency impact.
5. Create `.session/donna/023-8088de1701e0-notes.md` with evidence, Depmesh output, policy rationale, working-tree status, and exactly one of each:
   - `<!-- bec-port id=post-sync-021-8088de17-ignore-docker-readme status=changed -->` or `<!-- bec-port id=post-sync-021-8088de17-ignore-docker-readme status=no-change -->`
   - `<!-- focused-check id=post-sync-021-8088de17-ignore-docker-readme source-parse=passed test-parse=passed behavior=passed -->`
6. Then `{{ donna.lib.goto("verify_bec_port") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify the Source Decision

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
notes=".session/donna/023-8088de1701e0-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port id=post-sync-021-8088de17-ignore-docker-readme status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check id=post-sync-021-8088de17-ignore-docker-readme source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the Source Decision

```toml donna
id = "repair_bec_port"
kind = "donna.lib.request_action"
```

Verification failed. Review `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`, repair only this metadata decision or its notes after querying Depmesh, then `{{ donna.lib.goto("verify_bec_port") }}`; if unresolved, `{{ donna.lib.goto("blocked") }}`.

## Review the Source Decision

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report the evidence, changed/no-change decision, absence of dependency/test impact, checks, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or an explicit no-new-commit decision.
4. Only on that resume, record exactly one `<!-- source-review id=post-sync-021-8088de17-ignore-docker-readme commit=<40-hex-or-none> -->`.
5. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`; if corrections are requested, `{{ donna.lib.goto("repair_bec_port") }}`; if rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

## Reverify the Source Decision

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
notes=".session/donna/023-8088de1701e0-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review id=post-sync-021-8088de17-ignore-docker-readme commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Compare the Split Metadata

```toml donna
id = "compare_split_implementation"
kind = "donna.lib.request_action"
```

1. Use pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split` to read provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`.
2. Compare only `29a7478dc7deea846a2c1ff1abd25a881e6f87db:.gitignore`; it lacks `Docker_README.md` and adds `Rplots.pdf`.
3. Use pinned Git-object reads only; never inspect or write either external working tree and never access `man/**` or `NAMESPACE`.
4. Treat destination policy as authoritative and the stopped split as rejected evidence.
5. If comparison is sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Apply the Split Decision

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before any edit and apply only independently justified `.gitignore` policy; a no-change result is valid.
2. Preserve user changes and never stage, commit, amend, or touch generated files.
3. Update the notes with rationale and exactly one of each:
   - `<!-- split-stage id=post-sync-021-8088de17-ignore-docker-readme status=changed -->` or `<!-- split-stage id=post-sync-021-8088de17-ignore-docker-readme status=no-change -->`
   - `<!-- split-focused-check id=post-sync-021-8088de17-ignore-docker-readme source-parse=passed test-parse=passed behavior=passed -->`
4. Then `{{ donna.lib.goto("verify_split_adjustment") }}`; if blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify the Split Decision

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
notes=".session/donna/023-8088de1701e0-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-stage id=post-sync-021-8088de17-ignore-docker-readme status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check id=post-sync-021-8088de17-ignore-docker-readme source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do test -f "$path"; Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the Split Decision

```toml donna
id = "repair_split_adjustment"
kind = "donna.lib.request_action"
```

Review `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`, repair only this split decision or its notes after querying Depmesh, then `{{ donna.lib.goto("verify_split_adjustment") }}`; if unresolved, `{{ donna.lib.goto("blocked") }}`.

## Review the Split Decision

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report the comparison, changed/no-change decision, checks, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message.
3. Do not complete this action until the maintainer explicitly resumes and confirms a destination commit SHA or an explicit no-new-commit decision.
4. Only on that resume, record exactly one `<!-- split-review id=post-sync-021-8088de17-ignore-docker-readme commit=<40-hex-or-none> -->`.
5. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`; if corrections are requested, `{{ donna.lib.goto("repair_split_adjustment") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Reverify the Split Decision

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
notes=".session/donna/023-8088de1701e0-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-review id=post-sync-021-8088de17-ignore-docker-readme commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1
mapfile -t changed_r < <({ git diff --name-only -- R tests/testthat; git ls-files --others --exclude-standard -- R tests/testthat; } | sort -u | grep -E '^(R|tests/testthat)/.*[.]R$' || true)
for path in "${changed_r[@]}"; do Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null; done
Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Both review gates and reverification passed. Report the metadata decisions and maintainer-confirmed commit values to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved policy or verification issue, preserved working tree, and exact maintainer input required.
