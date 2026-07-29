# Establish Standalone Core Baseline

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Complete the migrated `proBatch` package as a standalone baseline, define the smallest public extension contract required by `proBatchBench`, restore Core-owned t-SNE and UMAP diagnostics from pinned original sources, and verify the editable and installed package surfaces without replaying the rejected split rewrite.

This workflow MUST NOT edit an external source or downstream repository, inspect or modify `man/`, stage files, create commits, reset the Donna session, or rename the package from `proBatch`.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, all specifications relevant to the artifacts in scope, and `{{ donna.lib.path("@/workflows/reports/bec-ecoli-core-migration-history.md") }}`.
2. Run `depmesh -p llm relations`, then query all relations for every project artifact that the work may change.
3. Run `donna -p llm status` and confirm there is no unrelated pending work. Do not reset the session.
4. Inspect the current worktree and preserve all unrelated user changes.
5. Verify these immutable local Git objects without reading dirty external working-tree files:
   - original proBatch plot addition `d2ae7da22d92a2c7de3d08d768c23f7582efc032`;
   - final original proBatch source `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92` in `/home/yuliya/repos/other/proBatch`;
   - rejected split comparator `29a7478dc7deea846a2c1ff1abd25a881e6f87db` in `/home/yuliya/repos/other/proBatch-core-split`;
   - downstream proBatchBench reference `60ace4572cacac1e5160ff4f147632427d7678c8` in `/home/yuliya/repos/cosybio/proBatchBench`.
6. Use `git show <sha>:<path>` or another object-level read for external evidence. Keep all external repositories read-only.
7. If the governed scope, archive, references, and user changes are understood, `{{ donna.lib.goto("specify_baseline_contract") }}`.
8. If an immutable input is unavailable or the scope conflicts with a current specification, `{{ donna.lib.goto("blocked") }}`.

## Specify Baseline Contract

```toml donna
id = "specify_baseline_contract"
kind = "donna.lib.request_action"
```

1. Add a compact behavior specification for the standalone Core baseline and update `specs/intro.md`, `depmesh.toml`, and affected `AGENTS.md` guidance in the same change.
2. Define `proBatch` as the independently installable package and `proBatchBench` as a one-way consumer. Core MUST NOT import Bench or provider-specific engines.
3. Assign Core ownership to general containers, transformations, correction primitives, missing-value policies, diagnostic plots, and PCA/t-SNE/UMAP over matrices and `ProBatchFeatures`.
4. Assign Bench ownership to provider integrations, benchmark orchestration, evaluation metrics, and downstream presentation that is not a general Core diagnostic.
5. Define the observable lineage, assay-name, sample-alignment, missing-value, and dependency contracts needed by the remaining fixes.
6. Inventory the pinned Bench checkout's actual `proBatch::` calls. Specify the smallest provider-neutral extension surface that makes a clean Bench implementation possible, including only capabilities proven necessary for registration lifecycle, matrix adaptation, or structured step output.
7. Treat the stopped split as evidence only. Do not copy its registry, matrix-adapter, structured-result, identifier, dependency, R-version, or versioning choices without an independently specified Core requirement.
8. Define t-SNE/UMAP compatibility, optional-backend policy, reproducibility, and missing-value defaults consistently with the Core plotting family.
9. If the contract is complete and testable, `{{ donna.lib.goto("verify_contract") }}`.
10. If a product decision cannot be made from the stated standalone/downstream goal, `{{ donna.lib.goto("blocked") }}`.

## Verify Contract

```toml donna
id = "verify_contract"
kind = "donna.lib.run_script"
save_stdout_to = "contract_stdout"
save_stderr_to = "contract_stderr"
goto_on_success = "review_contract"
goto_on_failure = "repair_contract"
timeout = 180
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

test -f specs/behavior/core_baseline.md
rg -q '^# ' specs/behavior/core_baseline.md
rg -q '^## Goal of the document$' specs/behavior/core_baseline.md
rg -q '^## Scope$' specs/behavior/core_baseline.md
rg -q 'specs/behavior/core_baseline.md' specs/intro.md

depmesh -p llm relations
depmesh -p llm dependencies @/specs/behavior/core_baseline.md
depmesh -p llm dependencies --relation governs @/specs/behavior/core_baseline.md
depmesh -p llm dependencies @/AGENTS.md
donna -p llm validate @/workflows/establish-standalone-core-baseline.donna.md
git diff --check -- AGENTS.md depmesh.toml specs workflows
```

## Repair Contract

```toml donna
id = "repair_contract"
kind = "donna.lib.request_action"
```

The contract verification failed.

Standard output:

```text
{{ donna.lib.task_variable("contract_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("contract_stderr") }}
```

Repair only the specification, index, guidance, relation, or workflow inconsistency shown above. Then `{{ donna.lib.goto("verify_contract") }}`.

If repair requires a new maintainer decision, `{{ donna.lib.goto("blocked") }}`.

## Review Contract

```toml donna
id = "review_contract"
kind = "donna.lib.request_action"
```

Pause for maintainer review before changing package behavior.

1. Inspect the exact in-scope unstaged diff and in-scope untracked files without modifying the Git index.
2. Confirm that the contract supports a standalone `proBatch`, a one-way `proBatchBench` import, and Core ownership of PCA/t-SNE/UMAP.
3. Confirm that every proposed extension API is justified by a general Core contract and actual pinned Bench use rather than by the rejected split alone.
4. Do not stage, commit, amend, or inspect generated manuals.
5. If the contract is accepted, `{{ donna.lib.goto("implement_package_hygiene") }}`.
6. If it needs revision, `{{ donna.lib.goto("repair_contract") }}`.
7. If the maintainer cannot accept either path, `{{ donna.lib.goto("blocked") }}`.

## Implement Package Hygiene

```toml donna
id = "implement_package_hygiene"
kind = "donna.lib.request_action"
```

1. Query Depmesh for every affected build, metadata, source, and test artifact.
2. Exclude `.pixi`, `pixi.toml`, and `pixi.lock` from R source-package builds.
3. Raise the declared testthat minimum to at least `3.2.0` and synchronize the direct Pixi constraint.
4. Replace the test-only `reshape2::dcast()` use with an independent base-R construction.
5. Replace the sole `plyr::arrange()` use with stable base ordering, then remove only the now-unneeded direct `plyr` and `r-reshape2` requirements. Accept that transitive packages may remain in the lock.
6. Consolidate duplicated `ProBatchFeatures` link tests into one authoritative home.
7. Find and close the graphics-device leak that creates `Rplots.pdf`; do not hide the artifact with an ignore rule.
8. Remove avoidable plot-test warnings and perform only evidence-backed dependency cleanup. Keep dependencies with proven runtime or test use.
9. Synchronize `DESCRIPTION`, `pixi.toml`, and `pixi.lock` without changing the package name or release version.
10. When the package and test hygiene changes are complete, `{{ donna.lib.goto("verify_package_hygiene") }}`.
11. If completion is blocked by dependency resolution or environment access, `{{ donna.lib.goto("blocked") }}`.

## Verify Package Hygiene

```toml donna
id = "verify_package_hygiene"
kind = "donna.lib.run_script"
save_stdout_to = "hygiene_stdout"
save_stderr_to = "hygiene_stderr"
goto_on_success = "implement_core_invariants"
goto_on_failure = "repair_package_hygiene"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

pixi run Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
while IFS= read -r path; do
    pixi run Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done < <(find R tests -type f -name '*.R' -print)

pixi run Rscript -e '
paths <- c(".pixi", "pixi.toml", "pixi.lock")
excluded <- vapply(paths, function(path) tools:::inRbuildignore(path, "."), logical(1))
if (!all(excluded)) stop("Build-visible development artifacts: ", paste(paths[!excluded], collapse = ", "))
'

if rg -n 'reshape2::dcast|plyr::arrange' R tests; then
    echo "Direct reshape2/plyr use remains" >&2
    exit 1
fi
test ! -e tests/testthat/Rplots.pdf
test ! -e Rplots.pdf
git diff --check -- . ':(exclude)man/**'
```

## Repair Package Hygiene

```toml donna
id = "repair_package_hygiene"
kind = "donna.lib.request_action"
```

The package-hygiene verification failed.

Standard output:

```text
{{ donna.lib.task_variable("hygiene_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("hygiene_stderr") }}
```

Diagnose the exact failure, preserve unrelated changes, and repair only the package-hygiene scope. Then `{{ donna.lib.goto("verify_package_hygiene") }}`.

If the failure requires network access or a maintainer dependency decision that is not authorized, `{{ donna.lib.goto("blocked") }}`.

## Implement Core Invariants

```toml donna
id = "implement_core_invariants"
kind = "donna.lib.request_action"
```

1. Query Depmesh for the container, missing-filter, correction, and affected test files.
2. Make existing-assay additions reject conflicting data, parent lineage, and virtual target reservations while keeping idempotent identical additions valid.
3. Make lineage edge identity include all stable origin fields and reject a second non-self parent for one result.
4. Apply one documented policy to explicit result-name collisions in `pb_filterNA()` and `pb_groupfilterNA()`; explicit names SHOULD error rather than be silently changed.
5. Delete the three shadowed correction definitions while retaining the later authoritative deprecated forwarding wrappers.
6. Add focused lineage, naming, compatibility-wrapper, and top-level symbol-uniqueness regressions.
7. Keep the change incremental in the current files; do not replace them with stopped-split versions.
8. When the invariants and tests are complete, `{{ donna.lib.goto("verify_core_invariants") }}`.
9. If the accepted contract cannot be implemented without a broader API change, `{{ donna.lib.goto("blocked") }}`.

## Verify Core Invariants

```toml donna
id = "verify_core_invariants"
kind = "donna.lib.run_script"
save_stdout_to = "core_stdout"
save_stderr_to = "core_stderr"
goto_on_success = "review_core_repairs"
goto_on_failure = "repair_core_invariants"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

for path in R/ProBatchFeatures.R R/pb_missing_filters.R R/correct_batch_effects.R; do
    pixi run Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done

pixi run Rscript -e 'devtools::test(filter = "ProBatchFeatures|pb_missing_helpers|correct_batch_effect", reporter = "summary")'
test ! -e tests/testthat/Rplots.pdf
git diff --check -- R tests/testthat
```

## Repair Core Invariants

```toml donna
id = "repair_core_invariants"
kind = "donna.lib.request_action"
```

The focused Core verification failed.

Standard output:

```text
{{ donna.lib.task_variable("core_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("core_stderr") }}
```

Diagnose and repair only the specified Core invariants or their focused tests, then `{{ donna.lib.goto("verify_core_invariants") }}`.

If the failure exposes an unresolved public contract, `{{ donna.lib.goto("blocked") }}`.

## Review Core Repairs

```toml donna
id = "review_core_repairs"
kind = "donna.lib.request_action"
```

1. Inspect the exact Core-invariant diff and new tests without modifying the index.
2. Confirm that the implementation follows the accepted specification, retains compatibility wrappers, and contains no provider or Bench implementation.
3. Do not stage or create a commit.
4. If the repairs are accepted, `{{ donna.lib.goto("design_companion_boundary") }}`.
5. If they need correction, `{{ donna.lib.goto("repair_core_invariants") }}`.
6. If the contract must change first, `{{ donna.lib.goto("repair_contract") }}`.

## Design Companion Boundary

```toml donna
id = "design_companion_boundary"
kind = "donna.lib.request_action"
```

1. Read the pinned proBatchBench source at `60ace4572cacac1e5160ff4f147632427d7678c8` without changing its repository.
2. Inventory every public Core call and classify it as already available, replaceable by a simpler existing Core API, or a missing provider-neutral extension contract.
3. Pay particular attention to provider registration metadata, aliases, availability checks, unregister lifecycle, matrix adaptation, and structured step output.
4. Compare the rejected split only after establishing the required behavior independently.
5. Update the Core behavior specification with the minimal accepted extension surface and explicit exclusions.
6. If no new Core API is needed, `{{ donna.lib.goto("verify_companion_boundary") }}`.
7. If a bounded provider-neutral API is required, `{{ donna.lib.goto("implement_companion_boundary") }}`.
8. If Core and Bench ownership cannot be reconciled from the stated goal, `{{ donna.lib.goto("blocked") }}`.

## Implement Companion Boundary

```toml donna
id = "implement_companion_boundary"
kind = "donna.lib.request_action"
```

1. Implement only the provider-neutral public contracts accepted in the baseline specification.
2. Extend the current registry incrementally; do not transplant the stopped split's files or broaden identifier policy incidentally.
3. Keep Core independent of proBatchBench and every provider engine.
4. Test the extension lifecycle with an in-package fake provider, including registration, availability, aliases, replay or lookup behavior where specified, cleanup, matrix orientation and identity, and structured artifacts where specified.
5. Document every public API in its Roxygen2 source. Do not generate documentation or edit `NAMESPACE`.
6. When the minimal boundary is implemented, `{{ donna.lib.goto("verify_companion_boundary") }}`.
7. If implementation requires a contract expansion, `{{ donna.lib.goto("repair_contract") }}`.

## Verify Companion Boundary

```toml donna
id = "verify_companion_boundary"
kind = "donna.lib.run_script"
save_stdout_to = "boundary_stdout"
save_stderr_to = "boundary_stderr"
goto_on_success = "restore_embeddings"
goto_on_failure = "repair_companion_boundary"
timeout = 600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

while IFS= read -r path; do
    pixi run Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done < <(find R tests/testthat -type f -name '*.R' -print)

if rg -n 'proBatchBench::|Imports:.*proBatchBench|Depends:.*proBatchBench' R DESCRIPTION; then
    echo "Core must not depend on proBatchBench" >&2
    exit 1
fi

pixi run Rscript -e 'devtools::test(filter = "extension.contract|registry|ProBatchFeatures", reporter = "summary")'
depmesh -p llm dependencies @/R/ProBatchFeatures.R
git diff --check -- R tests/testthat DESCRIPTION
```

## Repair Companion Boundary

```toml donna
id = "repair_companion_boundary"
kind = "donna.lib.request_action"
```

The companion-boundary verification failed.

Standard output:

```text
{{ donna.lib.task_variable("boundary_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("boundary_stderr") }}
```

Repair the provider-neutral Core contract or its fake-provider tests without editing proBatchBench or adding provider dependencies. Then `{{ donna.lib.goto("verify_companion_boundary") }}`.

If the failure shows that the accepted ownership contract is insufficient, `{{ donna.lib.goto("repair_contract") }}`.

## Restore Embeddings

```toml donna
id = "restore_embeddings"
kind = "donna.lib.request_action"
```

1. Query Depmesh for `R/proteome_wide_diagnostics.R`, `R/plot_helpers.R`, their focused tests, `DESCRIPTION`, Pixi files, and both main vignettes.
2. Read the original implementation from Git objects `d2ae7da22d92a2c7de3d08d768c23f7582efc032` and `e2bb18547c73f1c471fc1afcb3facbd8bea5fa92`.
3. Restore `plot_TSNE()` and `plot_UMAP()` generics plus their `default` and `ProBatchFeatures` methods. Restore only missing shared Plotly and multi-assay helpers; reuse current Core helpers wherever behavior already exists.
4. Preserve feature-row/sample-column orientation, annotation alignment, ordered multi-assay results, categorical colour and shape handling, t-SNE perplexity safeguards, static ggplot defaults, optional Plotly rendering, and explicit reproducibility controls.
5. Keep the original `Rtsne` and `umap` backends. Do not substitute another backend or copy `R/embedding_metrics.R` from proBatchBench.
6. Keep missing-value defaults consistent across PCA, t-SNE, and UMAP. Any coordinated default change must be specified, documented, and tested rather than silently copied from Bench.
7. Namespace backend calls, guard optional rendering, declare `Rtsne`, `umap`, and `plotly` in `Suggests`, and add direct Pixi dependencies so backend paths are exercised.
8. Add focused matrix and `ProBatchFeatures` tests for static and optional interactive output, validation, ordering, missing values, forwarded backend arguments, seed behavior, single/multi-assay results, and subplot behavior.
9. Make examples self-contained in Roxygen2 sources. Do not generate `man/` or edit `NAMESPACE`.
10. When the embedding surface is complete, `{{ donna.lib.goto("verify_embeddings") }}`.
11. If an optional backend is unavailable in the supported environment, `{{ donna.lib.goto("blocked") }}`.

## Verify Embeddings

```toml donna
id = "verify_embeddings"
kind = "donna.lib.run_script"
save_stdout_to = "embeddings_stdout"
save_stderr_to = "embeddings_stderr"
goto_on_success = "review_public_surface"
goto_on_failure = "repair_embeddings"
timeout = 900
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

pixi run Rscript -e 'parse(file = "R/proteome_wide_diagnostics.R")' >/dev/null
pixi run Rscript -e '
required <- c("Rtsne", "umap", "plotly")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing embedding test backends: ", paste(missing, collapse = ", "))
'

rg -q '^plot_TSNE <- function' R/proteome_wide_diagnostics.R
rg -q '^plot_TSNE.default <- function' R/proteome_wide_diagnostics.R
rg -q '^plot_TSNE.ProBatchFeatures <- function' R/proteome_wide_diagnostics.R
rg -q '^plot_UMAP <- function' R/proteome_wide_diagnostics.R
rg -q '^plot_UMAP.default <- function' R/proteome_wide_diagnostics.R
rg -q '^plot_UMAP.ProBatchFeatures <- function' R/proteome_wide_diagnostics.R

pixi run Rscript -e 'devtools::test(filter = "proteome_wide_diagnostics", reporter = "summary")'
test ! -e tests/testthat/Rplots.pdf
git diff --check -- R tests/testthat DESCRIPTION pixi.toml pixi.lock
```

## Repair Embeddings

```toml donna
id = "repair_embeddings"
kind = "donna.lib.request_action"
```

The embedding verification failed.

Standard output:

```text
{{ donna.lib.task_variable("embeddings_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("embeddings_stderr") }}
```

Repair only the Core-owned t-SNE/UMAP implementation, shared helpers, focused tests, or direct dependency metadata. Then `{{ donna.lib.goto("verify_embeddings") }}`.

If a failure requires a deliberate backend or compatibility change, `{{ donna.lib.goto("repair_contract") }}`.

## Review Public Surface

```toml donna
id = "review_public_surface"
kind = "donna.lib.request_action"
```

1. Inspect the exact extension and embedding diff without modifying the index.
2. Confirm that Core contains general APIs only, that the stopped split and Bench files were not copied wholesale, and that every new export has focused tests and Roxygen2 source.
3. Confirm that t-SNE and UMAP are Core-owned while benchmark metrics and provider implementations remain Bench-owned.
4. Do not stage, commit, generate documentation, or inspect `man/`.
5. If the public surface is accepted, `{{ donna.lib.goto("synchronize_project") }}`.
6. If embedding behavior needs repair, `{{ donna.lib.goto("repair_embeddings") }}`.
7. If the extension boundary needs repair, `{{ donna.lib.goto("repair_companion_boundary") }}`.
8. If the ownership contract needs revision, `{{ donna.lib.goto("repair_contract") }}`.

## Synchronize Project

```toml donna
id = "synchronize_project"
kind = "donna.lib.request_action"
```

1. Complete exact bidirectional Depmesh mappings after the final source/test layout, including batch-correction helpers, PBF input helpers, and proteome-wide diagnostics.
2. Update canonical Roxygen2 examples and descriptions without reading or generating `man/`.
3. Update `vignettes/proBatch.Rmd` and `vignettes/proBatchFeatures.Rmd` to use current correction policies and the Core t-SNE/UMAP APIs.
4. Update `NEWS` from the complete migration and baseline work, including breaking defaults, the extension contract, and restored embeddings.
5. Correct independently verified citation defects without coupling unrelated citation redesign.
6. Reconcile `DESCRIPTION`, Pixi metadata and lock, Roxygen imports, optional dependency guards, specifications, and tests.
7. Leave the release version decision to the maintainer; do not copy the stopped split's R-version or package-version changes.
8. Ensure no stale migration report is reintroduced and no generated documentation is accessed.
9. When editable project artifacts are synchronized, `{{ donna.lib.goto("verify_editable_package") }}`.
10. If a release or dependency choice requires maintainer direction, `{{ donna.lib.goto("blocked") }}`.

## Verify Editable Package

```toml donna
id = "verify_editable_package"
kind = "donna.lib.run_script"
save_stdout_to = "editable_stdout"
save_stderr_to = "editable_stderr"
goto_on_success = "review_release_sources"
goto_on_failure = "repair_editable_package"
timeout = 1800
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

pixi run Rscript -e 'read.dcf("DESCRIPTION")' >/dev/null
while IFS= read -r path; do
    pixi run Rscript -e 'parse(file = commandArgs(TRUE)[1])' "$path" >/dev/null
done < <(find R tests -type f -name '*.R' -print)

while IFS= read -r spec; do
    rg -q '^# ' "$spec"
    rg -q '^## Goal of the document$' "$spec"
    rg -q '^## Scope$' "$spec"
done < <(find specs -type f -name '*.md' -print)

depmesh -p llm relations
depmesh -p llm dependencies @/R/ProBatchFeatures.R
depmesh -p llm dependencies @/R/proteome_wide_diagnostics.R
depmesh -p llm dependencies @/tests/testthat/test-proteome_wide_diagnostics.R
depmesh -p llm dependencies --relation governs @/specs/behavior/core_baseline.md
depmesh -p llm dependencies @/workflows/establish-standalone-core-baseline.donna.md

donna -p llm validate --all
pixi run Rscript -e 'devtools::test(reporter = "summary")'

test ! -e Rplots.pdf
test ! -e tests/testthat/Rplots.pdf
if rg -n 'proBatchBench::|Imports:.*proBatchBench|Depends:.*proBatchBench' R DESCRIPTION; then
    echo "Core must not depend on proBatchBench" >&2
    exit 1
fi
git diff --check -- . ':(exclude)man/**'
```

## Repair Editable Package

```toml donna
id = "repair_editable_package"
kind = "donna.lib.request_action"
```

The complete agent-safe verification failed.

Standard output:

```text
{{ donna.lib.task_variable("editable_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("editable_stderr") }}
```

Classify the failure, query affected Depmesh relations, and repair only non-generated project artifacts in scope. Then `{{ donna.lib.goto("verify_editable_package") }}`.

If the failure is an environment limitation or requires inspecting `man/`, `{{ donna.lib.goto("blocked") }}`.

## Review Release Sources

```toml donna
id = "review_release_sources"
kind = "donna.lib.request_action"
```

1. Inspect the exact in-scope editable diff and untracked files without modifying the Git index.
2. Confirm that specifications, implementation, tests, Roxygen2 sources, vignettes, dependencies, lock data, and `NEWS` describe the same standalone baseline.
3. Confirm that no generated, cache, check, or external-repository artifact is included.
4. Confirm that the version choice is explicit and maintainer-owned.
5. Do not stage or create a commit.
6. If editable release sources are accepted, `{{ donna.lib.goto("maintainer_generation_and_checks") }}`.
7. If repairs are needed, `{{ donna.lib.goto("repair_editable_package") }}`.
8. If the maintainer defers the release decision, `{{ donna.lib.goto("blocked") }}`.

## Maintainer Generation and Checks

```toml donna
id = "maintainer_generation_and_checks"
kind = "donna.lib.request_action"
```

The agent MUST NOT perform this operation on the maintainer's behalf.

1. Ask the maintainer to regenerate Roxygen output and `NAMESPACE` manually.
2. Ask the maintainer to build both vignettes and run the package and Bioconductor checks that may inspect `man/`.
3. Ask the maintainer to verify default ComBat behavior and quantile normalization in an unrestricted environment.
4. Do not inspect, compare, lint, validate, or repair any `man/` result.
5. If the maintainer confirms generation and all required checks, `{{ donna.lib.goto("verify_import_surface") }}`.
6. If a non-`man/` source failure is reported, `{{ donna.lib.goto("repair_editable_package") }}`.
7. If a generated-manual issue or environment blocker remains, `{{ donna.lib.goto("blocked") }}`.

## Verify Import Surface

```toml donna
id = "verify_import_surface"
kind = "donna.lib.request_action"
```

1. Load the maintainer-generated package namespace without inspecting `man/`.
2. Verify that the accepted extension APIs and `plot_TSNE()`/`plot_UMAP()` generics are exported and that their `default` and `ProBatchFeatures` S3 methods are registered.
3. Exercise a minimal fake-provider lifecycle and matrix/`ProBatchFeatures` embedding smoke test through the public namespace.
4. Compare the installed Core surface with the pinned proBatchBench call inventory. Keep the downstream repository read-only.
5. Confirm that proBatchBench can remove its duplicated embedding implementation and import the Core functions, and record remaining Bench-only refactoring in the final handoff.
6. If the installed Core is a sufficient standalone import target, `{{ donna.lib.goto("finish") }}`.
7. If a Core extension contract is incomplete, `{{ donna.lib.goto("repair_companion_boundary") }}`.
8. If a Core embedding contract is incomplete, `{{ donna.lib.goto("repair_embeddings") }}`.
9. If only downstream Bench changes remain, `{{ donna.lib.goto("finish") }}`.
10. If ownership is still ambiguous, `{{ donna.lib.goto("blocked") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

The standalone Core baseline is complete. Report the accepted public and ownership contracts, Core fixes, restored t-SNE/UMAP behavior, verification evidence, maintainer confirmations, and the exact downstream-only changes still required in proBatchBench.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The standalone Core baseline cannot continue without maintainer input or an external-state change. Report completed operations, the exact blocker, affected artifacts, preserved user changes, and the decision or environment change required. Do not stage, commit, edit external repositories, or inspect generated manuals.
