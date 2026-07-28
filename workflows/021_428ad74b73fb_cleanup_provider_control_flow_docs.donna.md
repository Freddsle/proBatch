# Migrate Core Documentation and LOESS Cleanup

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Port only core-owned documentation and LOESS cleanup from mixed commit `428ad74b73fbc067de09d0a71c93395ae85eb51f`, excluding all provider and embedding hunks, with two maintainer review pauses.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-019-428ad74b-cleanup-provider-control-flow-docs` in `{{ donna.lib.path("@/workflows/bec-ecoli-core-migration-manifest.md") }}`.
2. Inspect commit `428ad74b73fbc067de09d0a71c93395ae85eb51f` against parent `6989542e4ece84f9fe210f3bf06b8254684e425e` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch`.
3. Non-generated paths are `R/M_ComBat.R`, `R/ProBatchFeatures.R`, `R/correct_batch_BERT_PLSDA.R`, `R/correct_batch_NormAE.R`, `R/correct_batch_RUVIIIC.R`, `R/correct_batch_effects_old.R`, `R/fit_non_linear.R`, `R/impute_PRONE.R`, `R/plot_NA_intensity.R`, `R/proBatch.R`, and `R/proteome_wide_diagnostics.R`.
4. Retain only core hunks: return documentation for `ProBatchFeatures`, the package topic, and legacy correction wrappers; environment-based LOESS warning capture; and `\donttest` for `plot_NA_intensity`.
5. Exclude at hunk level mComBat, PLSDA, NormAE, RUVIIIC, PRONE, TSNE, and UMAP cleanup and all provider dependencies.
6. Historical tests changed none. Focus on `test-ProBatchFeatures.R`, `test-correct_batch_effect.R`, `test-fit_non_linear.R`, and `test-plot_NA_intensity.R`; documentation checks remain source-only.
7. Split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` already has the NA-intensity `\donttest`; LOESS environment capture and some return docs are absent, and correction docs moved into consolidated `R/correct_batch_effects.R`.
8. Query Depmesh before editing. Never inspect external working trees, write externally, cherry-pick, or access `man/**` or `NAMESPACE`.
9. If ownership and relocated documentation targets are clear, `{{ donna.lib.goto("port_bec_behavior") }}`; otherwise `{{ donna.lib.goto("blocked") }}`.

## Port the BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

1. Compare current source comments and LOESS control flow before editing. Apply return docs to current Roxygen owners, not the removed `correct_batch_effects_old.R` path.
2. Evaluate LOESS warning capture in sequence after workflow 019; preserve equivalent or stronger error behavior. Keep the already-present NA-intensity example unchanged.
3. Run source parse plus focused class/correction/fit/plot tests as applicable. Keep every provider/embedding hunk out of core.
4. Use `apply_patch`, preserve user changes, and never cherry-pick, stage, commit, amend, draft a commit message, regenerate documentation, or touch generated files.
5. Create `.session/donna/021-core-docs-loess-cleanup-notes.md` with evidence, Depmesh/dependency and relocation decisions, checks, and exactly one of each:
   - `<!-- migration-id id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs -->`
   - `<!-- bec-port status=changed -->` or `<!-- bec-port status=no-change -->`
   - `<!-- focused-check source-parse=passed test-parse=passed behavior=passed -->`
6. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/021-core-docs-loess-cleanup-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs -->' "$notes")" -eq 1
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

Repair only retained core Roxygen/comments, LOESS behavior/tests, or notes after querying Depmesh. Preserve all provider exclusions. Then `{{ donna.lib.goto("verify_bec_port") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

This is the first hard maintainer review and commit pause.

1. Report retained/relocated/already-equivalent docs, LOESS result, exclusions, tests, and working-tree status.
2. Do not stage, commit, amend, or draft a commit message; the agent must not perform the maintainer's commit.
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
notes=".session/donna/021-core-docs-loess-cleanup-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs -->' "$notes")" -eq 1
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

1. Compare accepted behavior with commit `29a7478dc7deea846a2c1ff1abd25a881e6f87db` using only pinned Git-object commands in `/home/yuliya/repos/other/proBatch-core-split`.
2. Inspect `R/ProBatchFeatures.R`, `R/correct_batch_effects.R`, `R/fit_non_linear.R`, `R/plot_NA_intensity.R`, and `R/proBatch.R` plus focused tests. Confirm the plot example is equivalent and evaluate missing/relocated return docs and LOESS capture.
3. Confirm provider implementations, TSNE, and UMAP remain absent from core. Never inspect/write external working trees or generated files.
4. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}` if evidence conflicts.

## Apply the Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

1. Query Depmesh, then apply only independently justified core Roxygen/comment or LOESS hardening. Record each split target as adopted, equivalent, relocated, or excluded.
2. Preserve user changes; use `apply_patch`; never stage, commit, amend, cherry-pick, generate documentation, or touch generated files.
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
notes=".session/donna/021-core-docs-loess-cleanup-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs -->' "$notes")" -eq 1
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

Repair only retained core docs/LOESS changes, tests, or notes after querying Depmesh. Preserve all exclusions. Then `{{ donna.lib.goto("verify_split_adjustment") }}` or `{{ donna.lib.goto("blocked") }}`.

## Review the Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

This is the second hard maintainer review and commit pause.

1. Report split evidence, adopted/equivalent/relocated/excluded outcomes, checks, and working-tree status.
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
notes=".session/donna/021-core-docs-loess-cleanup-notes.md"
test -f "$notes"
test "$(grep -Fxc '<!-- migration-id id=post-sync-019-428ad74b-cleanup-provider-control-flow-docs -->' "$notes")" -eq 1
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

Report core documentation/LOESS outcomes, provider exclusions, and review markers from `.session/donna/021-core-docs-loess-cleanup-notes.md`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved ownership/documentation/API or verification decision, preserved working tree, and required maintainer input.
