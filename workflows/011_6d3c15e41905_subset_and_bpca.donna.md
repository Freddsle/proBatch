# Migrate Sample Subsetting and Review Bayesian PCA

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "inspect_source_commit"
```

Migrate core sample-subsetting behavior from source commit `6d3c15e419058d264b61ac43c27e95e0d635ca01` and make an explicit reviewed decision on optional Bayesian PCA.

## Inspect the Source Commit

```toml donna
id = "inspect_source_commit"
kind = "donna.lib.request_action"
```

1. Read manifest entry `post-sync-009-6d3c15e4-subset-and-bpca` and the ownership map.
2. Inspect pinned parent `0dbce8129444c9079a560fb9bcd60b328b1f054c`, commit `6d3c15e419058d264b61ac43c27e95e0d635ca01`, and diff in `/home/yuliya/repos/other/proBatch`.
3. Changed non-generated paths: `DESCRIPTION`, `R/ProBatchFeatures.R`, `R/plot_helpers.R`, `R/proteome_wide_diagnostics.R`, `tests/testthat/test-ProBatchFeatures.R`, `tests/testthat/test-plot_helpers.R`, and `tests/testthat/test-proteome_wide_diagnostics.R`.
4. Core hunks fix the four-index `ProBatchFeatures` `[` method, add `pb_subset_samples`, add shared/per-assay atomic argument handling, and add `.pb_compute_pca_embedding`/`.pb_format_pc_axis_label` plus BPCA parameters to `plot_PCA.default`.
5. Focused tests cover metadata subsetting, absent subset columns, three-index assay preservation, per-assay title vectors, and BPCA with missing data.
6. Dependency/API question: source adds optional `pcaMethods` to `Suggests`; decide explicitly whether BPCA belongs in lean core or is skipped. Do not infer acceptance from historical presence.
7. Exact split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`; targets: `R/ProBatchFeatures.R`, `R/plot_helpers.R`, `R/proteome_wide_diagnostics.R`, three focused tests, and `DESCRIPTION`. Subsetting exists; atomic mode/BPCA do not.
8. Skip rationale: do not duplicate accepted subsetting; BPCA and `pcaMethods` may be skipped only with a recorded core-policy decision.
9. Query Depmesh before edits. External reads are pinned Git objects only; never inspect/write external working trees or `man/**`/`NAMESPACE`. Never run `checkout`, `switch`, `reset`, `stash`, or `clean` in either external repository. These prohibitions apply to every operation in this workflow: preserve every pre-existing user change; never cherry-pick, stage, commit, amend, draft a commit message, or touch generated files.
10. Continue to `{{ donna.lib.goto("port_bec_behavior") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Port BEC Behavior

```toml donna
id = "port_bec_behavior"
kind = "donna.lib.request_action"
```

Compare each hunk independently; port missing core behavior/tests/Roxygen and minimum accepted metadata only. Preserve user changes; no cherry-pick, staging, commit, amend, commit-message drafting, or generated files. Create `.session/donna/011-6d3c15e4-subset-and-bpca-notes.md` with evidence, explicit BPCA dependency decision, Depmesh queries, checks, and exactly one:

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
notes=".session/donna/011-6d3c15e4-subset-and-bpca-notes.md"
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

Inspect `{{ donna.lib.task_variable("bec_verify_stdout") }}` and `{{ donna.lib.task_variable("bec_verify_stderr") }}`. Repair only this commit scope/notes after Depmesh queries. Then `{{ donna.lib.goto("verify_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Review BEC Port

```toml donna
id = "review_bec_port"
kind = "donna.lib.request_action"
```

Hard maintainer review and commit pause. Report subsetting, atomic-vector, BPCA/dependency decisions, checks, and status. The agent must not stage, commit, amend, or draft a commit message. Do not complete until explicit maintainer resume with a commit SHA or no-new-commit decision. Add exactly one `<!-- source-review commit=<40-hex-or-none> -->` marker. Then `{{ donna.lib.goto("reverify_bec_port") }}`, `{{ donna.lib.goto("repair_bec_port") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/011-6d3c15e4-subset-and-bpca-notes.md"
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

Use `git -C /home/yuliya/repos/other/proBatch-core-split show`, `diff`, `ls-tree`, or `grep` only against provenance `a04c9a29a3a4ba9f719d3b9c778616f3dd77903b:SPLIT_ATTEMPT.md` and exact comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db`. Reassess existing subsetting and missing atomic/BPCA behavior independently, including `DESCRIPTION`. Never write externally or inspect generated files. Then `{{ donna.lib.goto("apply_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

## Apply Split Adjustment

```toml donna
id = "apply_split_adjustment"
kind = "donna.lib.request_action"
```

After Depmesh queries, apply only independently justified hardening/tests/metadata consistent with the reviewed BPCA policy. Preserve user changes; no staging, commit, amend, or generated files. Add one concrete split-stage marker and `<!-- split-focused-check source-parse=passed test-parse=passed behavior=passed -->`. Then `{{ donna.lib.goto("verify_split_adjustment") }}`, or `{{ donna.lib.goto("blocked") }}`.

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
notes=".session/donna/011-6d3c15e4-subset-and-bpca-notes.md"
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
notes=".session/donna/011-6d3c15e4-subset-and-bpca-notes.md"
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

Report subsetting, atomic-vector, BPCA/dependency, and split decisions with both review markers.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

Report the unresolved core API/dependency decision and exact maintainer input required.
