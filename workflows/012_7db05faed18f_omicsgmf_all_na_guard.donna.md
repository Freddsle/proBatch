# Migrate Core colData Compatibility from omicsGMF Guard

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate only the core colData compatibility hunk from mixed source commit `7db05faed18f75572d92553ae668288c326bda2a`, excluding provider all-NA behavior.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-010-7db05fae-omicsgmf-all-na-guard` and the complete ownership map.
2. Inspect pinned parent `6d3c15e419058d264b61ac43c27e95e0d635ca01`, commit `7db05faed18f75572d92553ae668288c326bda2a`, and diff in `/home/yuliya/repos/other/proBatch`.
3. Changed non-generated paths: `R/ProBatchFeatures.R`, `R/impute_omicsGMF.R`, `tests/testthat/test-ProBatchFeatures.R`, `tests/testthat/test-impute_PRONE.R`, `tests/testthat/test-impute_omicsGMF.R`, and `tests/testthat/test-prone_normalization_steps.R`.
4. Retain core `.pb_harmonize_colData`: integer-like factor levels may equal matching numeric/integer values, while genuine mismatches still fail. Focused core tests cover both cases.
5. Exclude companion `.omicsgmf_fit_and_impute` all-NA removal/restoration and all omicsGMF/PRONE provider tests reached through `impute_with_omicsGMF`, `imputePRONE*`, and PRONE providers.
6. Dependency question: core uses existing container dependencies only; do not add `omicsGMF`, `SingleCellExperiment`, or `PRONE`. `SummarizedExperiment`/`S4Vectors` remain existing core infrastructure.
7. Exact split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; targets: `R/ProBatchFeatures.R`, `tests/testthat/test-ProBatchFeatures.R`, `R/registry.R`, `tests/testthat/test-registry.R`, and `tests/testthat/test-symbol-ownership.R`.
8. Skip rationale: provider algorithms and tests belong exclusively in the companion; retain only the independently useful core colData comparison branch/test.
9. Query Depmesh before edits. External reads are pinned Git objects only; never inspect/write external working trees or `man/**`/`NAMESPACE`. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
10. Continue to `{{ donna.lib.goto("port_bec_behavior") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

Compare first; port only `.pb_harmonize_colData`, Roxygen if applicable, and its two core tests. Preserve user changes; no cherry-pick, staging, commit, amend, commit-message drafting, provider code, or generated files. Create `.session/donna/012-7db05fae-omicsgmf-all-na-guard-notes.md` with hunk exclusions, dependency decision, Depmesh queries, checks, and exactly one:

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
notes=".session/donna/012-7db05fae-omicsgmf-all-na-guard-notes.md"
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

Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only core colData behavior/tests/notes after Depmesh queries; never import providers. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. Report core branch, provider exclusions, dependencies, checks, and status. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_bec_port") }}`, `{{ donna.lib.goto("repair_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/012-7db05fae-omicsgmf-all-na-guard-notes.md"
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

Use `git -C /home/yuliya/repos/other/proBatch-core-split show`, `diff`, `ls-tree`, or `grep` only against provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and exact comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`. Compare `.pb_harmonize_colData`, its focused test, and registry/negative ownership tests. Provider source absence is intentional. Never write externally or inspect generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

After Depmesh queries, apply only independently justified core compatibility/tests or negative ownership hardening. Preserve user changes; no staging, commit, amend, provider code, or generated files. Add one concrete split-stage marker and `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/012-7db05fae-omicsgmf-all-na-guard-notes.md"
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

Inspect `{{ donna.lib.task_variable("split_verify_stdout") }}` and `{{ donna.lib.task_variable("split_verify_stderr") }}`. Repair only this core adjustment after Depmesh queries. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review Split Adjustment

```toml donna
id = "review_split_adjustment"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- split-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_split_adjustment") }}`, `{{ donna.lib.goto("repair_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/012-7db05fae-omicsgmf-all-na-guard-notes.md"
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

Report core colData, provider-exclusion, and split decisions with both review markers.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved compatibility, ownership, dependency, or verification decision.
