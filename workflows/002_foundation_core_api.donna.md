# Migrate the Foundation Core API

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate the core-owned API present at synchronizing merge `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab`, then compare it with the rejected stopped-split implementation under two maintainer-controlled review pauses.

## Inspect the Foundation Snapshot

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read the `foundation-core-api` entry in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}` and the complete ownership map at `/home/yuliya/repos/other/proBatch/archive/stopped-probatch-split-2026-07-27/API_FUNCTIONS.md`.
2. Inspect the immutable source tree and diff at merge `5cd1d2f373ae3f515bad4702d6524fe44e1ed8ab` against baseline `ba6ee246eace090e71baa7aba302ca64e76ddb32` in `/home/yuliya/repos/other/proBatch`. Use only `git -C ... show`, `diff`, `ls-tree`, or `grep` against full commit IDs, or a read-only archive under `/tmp`.
3. Never inspect either external working tree and never run checkout, switch, reset, stash, clean, or any command that writes to an external repository. Never read, enumerate, compare, or transfer `man/**` or generated `NAMESPACE`.
4. Recheck all 104 foundation exports and their Roxygen sources. Resolve the effective definitions of `correct_batch_effects_df`, `correct_batch_effects_dm`, and `correct_with_removeBatchEffect_dm` from the later-collated `R/correct_batch_effects_old.R`; do not choose authority by filename recency.
5. Exclude the 26 companion exports, provider implementations/registrations, benchmark code, and companion-only dependencies at hunk level. Retain independent core PVCA and diagnostic hunks from mixed files.
6. Expected destination source families are `DESCRIPTION`, `NEWS`, `R/ProBatchFeatures.R`, `R/auxiliary.R`, `R/CV_calculation.R`, `R/colors_for_annotation.R`, `R/correct_batch_effects.R`, `R/correlation-based_diagnostics.R`, `R/date_conversion.R`, `R/design_diagnostics.R`, `R/feature_level_diagnostics.R`, `R/fit_non_linear.R`, `R/handle_missing_values.R`, `R/initial_assessment.R`, `R/metadata_diagnostics.R`, `R/normalize.R`, `R/pb_missing_filters.R`, `R/pbf_input_helpers.R`, `R/plot_helpers.R`, `R/plot_missing.R`, `R/plot_split_violinplot.R`, `R/proBatch.R`, `R/proteome_wide_diagnostics.R`, `R/transform_raw.R`, `R/utility_funcs.R`, and `R/zzz_helpers.R`. Add a destination file only when its retained core behavior warrants it.
7. Focused tests are the core-owned cases listed in the manifest, especially container/lineage, correction, missing-data, design/metadata, normalization, PVCA, and diagnostic tests.
8. Query all Depmesh relations for every destination artifact before editing it.
9. If ownership and the effective source behavior are clear, `{{ donna.lib.goto("port_bec_behavior") }}`.
10. If an ownership, load-order, or dependency decision cannot be resolved from the pinned evidence, `{{ donna.lib.goto("blocked") }}`.

## Port the Foundation Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare every retained core behavior with the current destination and detect already-equivalent code before editing.
2. Port only missing core behavior, its editable Roxygen comments under `R/`, minimum package metadata, and focused tests. Never copy a source file wholesale when it contains companion-owned hunks.
3. Preserve unrelated user changes. Do not cherry-pick, stage, commit, amend, draft a commit message, or modify generated files.
4. Use `apply_patch` for project edits. Leave documentation generation and `NAMESPACE` regeneration to the maintainer.
5. Run the narrowest agent-safe focused checks that do not read or validate `man/**` or `NAMESPACE`.
6. Create or update `.session/donna/002-foundation-core-api-notes.md` with evidence, dependency queries, changed/no-change decisions, focused check commands/results, and exactly one marker of each form:
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
7. When the source port and notes are ready, `{{ donna.lib.goto("verify_bec_port") }}`.
8. If safe migration is blocked, `{{ donna.lib.goto("blocked") }}`.

## Verify the Foundation Port

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

notes=".session/donna/002-foundation-core-api-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- bec-port status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1

mapfile -t changed_r < <(
    {
        git diff --name-only -- R tests/testthat
        git ls-files --others --exclude-standard -- R tests/testthat
    } |
        sort -u |
        grep -E '^(R|tests/testthat)/.*[.]R$' || true
)

for path in "${changed_r[@]}"; do
    test -f "$path"
    Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done

Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair the Foundation Port

```toml donna
id = "repair_bec_port"
kind = "donna.lib.request_action"
```

Foundation source-port verification failed.

Standard output:

```text
{{ donna.lib.task_variable("bec_verify_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("bec_verify_stderr") }}
```

1. Diagnose and repair only the foundation port, focused tests, metadata, or session notes. Query Depmesh before any additional governed edit.
2. Preserve companion exclusions and the generated-file prohibition.
3. When repaired, `{{ donna.lib.goto("verify_bec_port") }}`.
4. If repair requires an unresolved API or ownership decision, `{{ donna.lib.goto("blocked") }}`.

## Review the Foundation Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first maintainer-controlled review and commit pause.

1. Report the source evidence, destination changes or no-change decision, dependency impact, focused verification, and working-tree status to the maintainer.
2. Give the maintainer time to inspect, edit, stage, and manually commit if desired.
3. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
4. Do not complete this action request until the maintainer explicitly asks to resume and confirms either the resulting destination commit SHA or an explicit no-new-commit decision.
5. On explicit resume, record exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker in the session notes.
6. If accepted, `{{ donna.lib.goto("reverify_bec_port") }}`.
7. If corrections are requested, `{{ donna.lib.goto("repair_bec_port") }}`.
8. If the port is rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

## Reverify the Foundation Port

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

notes=".session/donna/002-foundation-core-api-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- source-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1

mapfile -t destination_r < <(find R tests/testthat -type f -name '*.R' -print | sort)
for path in "${destination_r[@]}"; do
    Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done

Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Compare the Stopped Split

```toml donna
id = "compare_split_implementation"
kind = "donna.lib.request_action"
```

1. Read the pinned provenance note with `git -C /home/yuliya/repos/other/proBatch-core-split show a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md`.
2. Compare the completed foundation behavior with stopped split commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db`, using pinned Git object commands only.
3. Inspect the split targets named in the manifest: container/lineage, registry, identifier, matrix-adapter, structured-result, correction, missing-data, design/metadata, diagnostics, tests, and minimum dependency metadata.
4. Treat the split as rejected evidence. Reassess each hardening change; do not copy it wholesale.
5. Do not inspect generated `man/**` or `NAMESPACE`, and do not write to the split repository.
6. If comparison evidence is sufficient, `{{ donna.lib.goto("apply_split_adjustment") }}`.
7. If a core API or ownership decision remains unresolved, `{{ donna.lib.goto("blocked") }}`.

## Apply Justified Split Hardening

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query all Depmesh relations for every destination artifact before editing.
2. Apply only independently justified core hardening, focused tests, and Roxygen source changes. Record why every reviewed split target is adopted, already equivalent, deferred to the residual workflow, or excluded.
3. Preserve unrelated changes; never stage, commit, amend, or touch generated files.
4. Run the narrowest agent-safe focused checks.
5. Update `.session/donna/002-foundation-core-api-notes.md` with exactly one marker of each form:
   - `<!-- split-stage status=changed -->` or `<!-- split-stage status=no-change -->`
   - `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`
6. When ready, `{{ donna.lib.goto("verify_split_adjustment") }}`.
7. If the split comparison exposes a blocker, `{{ donna.lib.goto("blocked") }}`.

## Verify Split Hardening

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

notes=".session/donna/002-foundation-core-api-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-stage status=(changed|no-change) -->$' "$notes")" -eq 1
test "$(grep -Ec '^<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->$' "$notes")" -eq 1

mapfile -t changed_r < <(
    {
        git diff --name-only -- R tests/testthat
        git ls-files --others --exclude-standard -- R tests/testthat
    } |
        sort -u |
        grep -E '^(R|tests/testthat)/.*[.]R$' || true
)

for path in "${changed_r[@]}"; do
    test -f "$path"
    Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done

Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Repair Split Hardening

```toml donna
id = "repair_split_adjustment"
kind = "donna.lib.request_action"
```

Split-stage verification failed.

Standard output:

```text
{{ donna.lib.task_variable("split_verify_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("split_verify_stderr") }}
```

1. Repair only the justified split-stage changes, focused tests, or session notes. Query Depmesh before additional edits.
2. Do not reopen accepted source-port behavior without reporting the conflict.
3. When repaired, `{{ donna.lib.goto("verify_split_adjustment") }}`.
4. If repair needs an unresolved ownership/API decision, `{{ donna.lib.goto("blocked") }}`.

## Review Split Hardening

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second maintainer-controlled review and commit pause.

1. Report the split comparison, adopted or rejected hardening, destination changes or no-change decision, focused verification, and working-tree status.
2. Give the maintainer time to inspect, edit, stage, and manually commit if desired.
3. Do not stage, commit, or amend; the maintainer performs the commit. Inspect the exact in-scope unstaged diff and any in-scope untracked files without modifying the index. If the change set is nonempty, validate it against `@/specs/general/commits.md` and provide the complete copy-paste commit message in a fenced `text` block in the chat handoff. If one commit cannot accurately describe the change set, provide one validated message per coherent commit and identify each intended file set. If the change set is empty, state that no commit message is needed.
4. Do not complete this action request until the maintainer explicitly asks to resume and confirms either the resulting destination commit SHA or an explicit no-new-commit decision.
5. On explicit resume, record exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker in the session notes.
6. If accepted, `{{ donna.lib.goto("reverify_split_adjustment") }}`.
7. If corrections are requested, `{{ donna.lib.goto("repair_split_adjustment") }}`.
8. If rejected or blocked, `{{ donna.lib.goto("blocked") }}`.

## Reverify Split Hardening

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

notes=".session/donna/002-foundation-core-api-notes.md"
test -f "$notes"
test "$(grep -Ec '^<!-- split-review commit=(none|[0-9a-f]{40}) -->$' "$notes")" -eq 1

mapfile -t destination_r < <(find R tests/testthat -type f -name '*.R' -print | sort)
for path in "${destination_r[@]}"; do
    Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done

Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
git diff --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
git diff --cached --check -- . ':(exclude)man/**' ':(exclude)NAMESPACE'
```

## Finish Foundation Migration

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The foundation workflow has passed both post-review verification gates. Report the exact source/split change decisions and the developer-confirmed source and split commit SHAs—or `none`—from `.session/donna/002-foundation-core-api-notes.md` to the parent workflow.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The foundation cannot continue safely. Report the unresolved API, ownership, dependency, or verification decision; the preserved working tree; and the exact maintainer input required to the parent workflow.
