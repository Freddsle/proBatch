# Review Version and variancePartition Metadata

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Audit metadata commit `72d11d1f7f92a10cb9a5244e68f9c19e830197b9`, retain core release ownership, exclude the companion dependency, and compare pinned split evidence under two maintainer-controlled pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-020-72d11d1f-version-variancepartition` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `72d11d1f7f92a10cb9a5244e68f9c19e830197b9` against parent `428ad74b73fbc067de09d0a71c93395ae85eb51f` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. The sole non-generated path is `DESCRIPTION`. There are no R functions or focused runtime tests.
4. The hunks change version `1.99.5` to `2.1.0` and relax `variancePartition` from `>= 1.40.1` to `>= 1.36.2`.
5. Current destination metadata already reports `2.1.0`; determine whether this is independently authoritative and record no change unless the maintainer makes a release decision. `variancePartition` supports Bench-owned variance-partition APIs and must remain absent from core.
6. Split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` retains version `1.99.5` and has no `variancePartition`; neither hunk should be copied mechanically.
7. Query all Depmesh relations for `DESCRIPTION` before editing. Never inspect external working trees, write externally, cherry-pick, stage, commit, amend, or access `man/**` or `NAMESPACE`.
8. If the version/dependency disposition is clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Read current `DESCRIPTION`, validate it as DCF, and record whether version `2.1.0` is already equivalent. Do not alter the version without explicit release authority.
2. Confirm `variancePartition` remains absent from core Imports/Suggests and record its companion ownership. Make a metadata change only if independently required by current core.
3. Use `apply_patch`, preserve user changes, and never cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
4. Create `.session/donna/022-version-variancepartition-notes.md` with source evidence, Depmesh/release/dependency decisions, DCF checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-020-72d11d1f-version-variancepartition -->`
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
notes=".session/donna/022-version-variancepartition-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-020-72d11d1f-version-variancepartition -->' "$notes")" -eq 1
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

Repair only `DESCRIPTION` or the notes after querying Depmesh. Do not add the companion dependency or make an unauthorized release change. Then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report version equivalence/release authority, dependency exclusion, DCF verification, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
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
notes=".session/donna/022-version-variancepartition-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-020-72d11d1f-version-variancepartition -->' "$notes")" -eq 1
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

1. Compare `DESCRIPTION` at commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split`.
2. Record that the comparator's `1.99.5` is stale evidence rather than release authority and that its omission of `variancePartition` agrees with core ownership.
3. Never inspect/write external working trees, copy metadata wholesale, or access generated files.
4. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if release/dependency evidence conflicts.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh before any edit. Normally record no split-stage change; adjust core metadata only with independent current requirements and explicit release authority.
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
notes=".session/donna/022-version-variancepartition-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-020-72d11d1f-version-variancepartition -->' "$notes")" -eq 1
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

Repair only authorized metadata or notes after querying Depmesh; keep `variancePartition` absent from core. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split metadata evidence, release/dependency decisions, DCF verification, and working-tree status.
2. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
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
notes=".session/donna/022-version-variancepartition-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-020-72d11d1f-version-variancepartition -->' "$notes")" -eq 1
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

Report version/release and dependency ownership decisions plus review markers from `.session/donna/022-version-variancepartition-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved release/dependency or verification decision, preserved working tree, and required maintainer input.
