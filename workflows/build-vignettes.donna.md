# Build Package Vignettes

```toml donna
id = "primary"
kind = "donna.lib.workflow"
start_operation_id = "preflight"
```

Build every declared package vignette from editable package sources in isolated clean R processes, repair in-scope failures, and prove that generated rendering artifacts do not enter the source tree.

This workflow MUST NOT access `man/`, generate documentation, install or update packages or the Pixi environment, use the network, run a package-wide build or check, stage files, create commits, reset the Donna session, or modify generated `NAMESPACE`.

## Preflight

```toml donna
id = "preflight"
kind = "donna.lib.request_action"
```

1. Read `{{ donna.lib.path("@/AGENTS.md") }}`, `{{ donna.lib.path("@/specs/intro.md") }}`, `{{ donna.lib.path("@/specs/meta/general.md") }}`, `{{ donna.lib.path("@/specs/behavior/files_relations.md") }}`, and `{{ donna.lib.path("@/specs/general/workflows.md") }}` completely.
2. Read the built-in Donna workflow and Depmesh usage instructions, run `depmesh -p llm relations`, and query all configured relations separately for this workflow, its catalog, `{{ donna.lib.path("@/AGENTS.md") }}`, and every vignette, supporting asset, package source, data artifact, or metadata file already known to be in the initiating scope.
3. Inspect the non-generated worktree without modifying the Git index. Preserve unrelated user changes and distinguish the initiating scope from pre-existing failures.
4. Confirm this action belongs to the current Donna task and no unrelated action request is being displaced. Do not reset the session.
5. Confirm that this workflow validated and rendered correctly in view mode before execution. Run those read-only checks now if they were not already completed while no Donna operation was executing.
6. Confirm that vignette rendering is authorized but package installation, dependency installation, network access, documentation generation, package-wide building or checking, staging, and commits are not.
7. Do not read, enumerate, compare, hash, validate, or otherwise inspect any path under `man/`. Do not modify or regenerate `NAMESPACE`.
8. If the governed scope, current diff, and rendering boundary are understood, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
9. If specifications conflict or required context is unavailable, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Workflow Artifact

```toml donna
id = "verify_workflow_artifact"
kind = "donna.lib.run_script"
save_stdout_to = "workflow_stdout"
save_stderr_to = "workflow_stderr"
goto_on_success = "verify_vignette_environment"
goto_on_failure = "repair_workflow_artifact"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

workflow_path="workflows/build-vignettes.donna.md"
workflow_id="@/${workflow_path}"

test -f "$workflow_path"

artifact_entry_count="$(
    rg -Fxc -- '- Artifact: `./workflows/build-vignettes.donna.md`' \
        specs/general/workflows.md
)"
if (( artifact_entry_count != 1 )); then
    printf 'Workflow catalog must contain exactly one build-vignettes artifact entry; found %d.\n' \
        "$artifact_entry_count" >&2
    exit 1
fi

catalog_state="$(
    awk '
        /^### Build package vignettes$/ {
            in_section = 1
            next
        }
        in_section && /^### / {
            exit
        }
        in_section && /^- State: / {
            sub(/^- State: /, "")
            print
            exit
        }
    ' specs/general/workflows.md
)"
case "$catalog_state" in
    "implementation in progress"|"implemented") ;;
    *)
        printf 'Unexpected build-vignettes catalog state: %s\n' \
            "$catalog_state" >&2
        exit 1
        ;;
esac

if ! rg -Fq -- '@/workflows/build-vignettes.donna.md' AGENTS.md; then
    printf 'AGENTS.md must route vignette-build triggers to this workflow.\n' >&2
    exit 1
fi

relation_artifacts() {
    local relation="$1"
    local artifact="$2"
    local output
    output="$(depmesh -p llm dependencies --relation "$relation" "$artifact")"
    if printf '%s\n' "$output" | rg -Fqx -- '## warnings'; then
        printf 'Depmesh returned warnings for %s on %s:\n%s\n' \
            "$relation" "$artifact" "$output" >&2
        return 1
    fi
    printf '%s\n' "$output" |
        sed -n 's/^- \(@\/[^[:space:]]\+\)$/\1/p' |
        LC_ALL=C sort
}

actual_governance="$(relation_artifacts governed_by "$workflow_id")"
expected_governance="$(
    printf '%s\n' \
        @/specs/behavior/files_relations.md \
        @/specs/general/workflows.md |
        LC_ALL=C sort
)"
if [[ "$actual_governance" != "$expected_governance" ]]; then
    printf 'Unexpected workflow governance.\nExpected:\n%s\nActual:\n%s\n' \
        "$expected_governance" "$actual_governance" >&2
    exit 1
fi

for governing_spec in \
    @/specs/behavior/files_relations.md \
    @/specs/general/workflows.md
do
    if ! relation_artifacts governs "$governing_spec" |
        rg -Fqx -- "$workflow_id"; then
        printf 'Reverse governance does not contain %s for %s.\n' \
            "$workflow_id" "$governing_spec" >&2
        exit 1
    fi
done

if ! rg -Fqx -- '^workflows$' .Rbuildignore; then
    printf 'Project workflows must remain excluded from the source package.\n' >&2
    exit 1
fi

safe_diff_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks diff --check -- "${safe_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_diff_paths[@]}"

for text_path in \
    "$workflow_path" \
    AGENTS.md \
    specs/general/workflows.md
do
    if rg -n '[[:blank:]]+$' "$text_path"; then
        printf '%s contains trailing whitespace.\n' "$text_path" >&2
        exit 1
    fi
    if LC_ALL=C rg -n $'\r$' "$text_path"; then
        printf '%s contains CRLF line endings.\n' "$text_path" >&2
        exit 1
    fi
    if [[ -n "$(tail -c 1 -- "$text_path")" ]]; then
        printf '%s must end with a newline.\n' "$text_path" >&2
        exit 1
    fi
done

printf 'Workflow source, lifecycle state, agent routing, package exclusion, governance, and non-generated diff hygiene are synchronized.\n'
```

## Repair Workflow Artifact

```toml donna
id = "repair_workflow_artifact"
kind = "donna.lib.request_action"
```

Workflow-artifact verification failed.

Standard output:

```text
{{ donna.lib.task_variable("workflow_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("workflow_stderr") }}
```

1. Query all configured relations separately for every artifact to be edited.
2. Repair only this workflow, its catalog entry, affected agent guidance, its scalable governance relation, or an in-scope formatting discrepancy.
3. Preserve unrelated changes and the Git index. Do not add a workflow-specific Depmesh rule when the scalable workflow rule applies.
4. Keep the catalog state at `implementation in progress` until the first representative execution finishes.
5. Validate every control-flow change while Donna is awaiting this action.
6. After repair, `{{ donna.lib.goto("verify_workflow_artifact") }}`.
7. If repair requires a new contract or developer decision, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Vignette Environment

```toml donna
id = "verify_vignette_environment"
kind = "donna.lib.run_script"
save_stdout_to = "environment_stdout"
save_stderr_to = "environment_stderr"
goto_on_success = "review_vignette_scope"
goto_on_failure = "diagnose_environment"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

test -f DESCRIPTION
test -f NAMESPACE
test -f pixi.toml
test -f pixi.lock
test -d R
test -d vignettes
command -v pixi >/dev/null 2>&1

mapfile -t vignette_sources < <(
    find vignettes -maxdepth 1 -type f -name '*.Rmd' -print |
        LC_ALL=C sort
)
vignette_count="$(
    printf '%s\n' "${vignette_sources[@]}" |
        wc -l |
        tr -d ' '
)"
if (( vignette_count == 0 )); then
    printf 'No vignette sources were found.\n' >&2
    exit 1
fi

generated_paths="$(
    find vignettes -mindepth 1 \
        \( \
            -type f \
            \( -name '*.html' -o -name '*.R' -o -name '.Rhistory' \) \
            -o -type d \
            \( -name '*_cache' -o -name '*_files' \) \
        \) -print
)"
if [[ -n "$generated_paths" ]]; then
    printf 'Generated vignette artifacts already exist in the source tree:\n%s\n' \
        "$generated_paths" >&2
    exit 1
fi

status_before="$(mktemp /tmp/probatch-vignette-status-before.XXXXXXXX)"
status_after="$(mktemp /tmp/probatch-vignette-status-after.XXXXXXXX)"
environment_script="$(mktemp /tmp/probatch-vignette-environment.XXXXXXXX.R)"
cleanup() {
    rm -f -- "$status_before" "$status_after" "$environment_script"
}
trap cleanup EXIT

safe_status_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
    "${safe_status_paths[@]}" > "$status_before"
pixi_before="$(sha256sum pixi.toml pixi.lock)"

cat > "$environment_script" <<'RSCRIPT'
description <- read.dcf("DESCRIPTION")

parse_dependencies <- function(field) {
    if (!field %in% colnames(description) || is.na(description[1L, field])) {
        return(character())
    }
    entries <- trimws(unlist(strsplit(description[1L, field], ",", fixed = TRUE)))
    trimws(sub("[[:space:]]*\\(.*\\)[[:space:]]*$", "", entries))
}

declared <- unique(unlist(lapply(
    c("Imports", "Suggests", "VignetteBuilder"),
    parse_dependencies
)))
required_declared <- c("BiocStyle", "knitr", "rmarkdown")
undeclared <- setdiff(required_declared, declared)
if (length(undeclared)) {
    stop("Missing declared vignette dependencies: ", paste(undeclared, collapse = ", "))
}

required_runtime <- unique(c(declared, "pkgload", "BiocParallel"))
missing_runtime <- required_runtime[
    !vapply(required_runtime, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_runtime)) {
    stop("Unavailable declared/runtime packages: ", paste(missing_runtime, collapse = ", "))
}
if (!rmarkdown::pandoc_available()) {
    stop("Pandoc is unavailable")
}

sources <- sort(list.files(
    "vignettes",
    pattern = "[.]Rmd$",
    full.names = TRUE
))
if (!length(sources)) {
    stop("No vignette sources found")
}
for (source in sources) {
    lines <- readLines(source, warn = FALSE, encoding = "UTF-8")
    if (!any(grepl("%\\\\VignetteIndexEntry\\{", lines))) {
        stop("Missing VignetteIndexEntry in ", source)
    }
    if (!any(grepl("%\\\\VignetteEngine\\{knitr::rmarkdown\\}", lines))) {
        stop("Unexpected or missing vignette engine in ", source)
    }
}

pkgload::load_all(
    ".",
    compile = FALSE,
    attach = TRUE,
    export_all = FALSE,
    export_imports = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE,
    quiet = TRUE,
    warn_conflicts = TRUE,
    debug = FALSE
)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

cat(
    sprintf(
        "Validated %d vignette sources, declared dependencies, Pandoc %s, and source loading.\n",
        length(sources),
        as.character(rmarkdown::pandoc_version())
    )
)
RSCRIPT

pixi run --as-is Rscript --vanilla "$environment_script"

git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
    "${safe_status_paths[@]}" > "$status_after"
cmp -s -- "$status_before" "$status_after" || {
    printf 'Environment verification changed the non-generated worktree.\n' >&2
    exit 1
}
pixi_after="$(sha256sum pixi.toml pixi.lock)"
if [[ "$pixi_before" != "$pixi_after" ]]; then
    printf 'Environment verification changed Pixi metadata.\n' >&2
    exit 1
fi
```

## Diagnose Environment

```toml donna
id = "diagnose_environment"
kind = "donna.lib.request_action"
```

Vignette-environment verification failed.

Standard output:

```text
{{ donna.lib.task_variable("environment_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("environment_stderr") }}
```

1. Determine whether the failure is an in-scope metadata, source, asset, workflow, or current-environment issue.
2. Query all configured relations separately before editing any governed artifact.
3. Do not install or update packages, modify the Pixi environment, use the network, inspect generated manuals, or modify `NAMESPACE`.
4. If an authorized source or workflow repair can make the environment valid, make that focused repair and `{{ donna.lib.goto("verify_vignette_environment") }}`.
5. If the environment cannot satisfy the declared dependencies without installation or another maintainer action, `{{ donna.lib.goto("record_blocker") }}`.

## Review Vignette Scope

```toml donna
id = "review_vignette_scope"
kind = "donna.lib.request_action"
```

The workflow artifact and local vignette environment passed their deterministic gates.

1. Inspect every declared vignette source and its supporting assets without accessing generated manuals.
2. Inspect the package APIs, data, metadata, and declared dependencies used by executable chunks. Query Depmesh separately for every artifact that may need repair.
3. Confirm installation examples and any network-using examples are non-evaluated.
4. Confirm the representative build may use the current editable sources loaded by `pkgload`, copied into an isolated temporary tree, with each vignette rendered by `rmarkdown` in its own `pixi run --as-is Rscript --vanilla` process.
5. If the scope and offline build plan are sound, `{{ donna.lib.goto("run_vignette_build") }}`.
6. If a source, asset, metadata, or workflow repair is needed first, make the focused repair and `{{ donna.lib.goto("verify_vignette_environment") }}`.
7. If a product decision or unavailable dependency blocks the build, `{{ donna.lib.goto("record_blocker") }}`.

## Run Vignette Build

```toml donna
id = "run_vignette_build"
kind = "donna.lib.run_script"
save_stdout_to = "build_stdout"
save_stderr_to = "build_stderr"
goto_on_success = "review_build_results"
goto_on_failure = "diagnose_build_failure"
timeout = 3600
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

project_root="$(pwd -P)"
build_root="$(mktemp -d /tmp/probatch-build-vignettes.XXXXXXXX)"
stage_root="$build_root/package"
render_script="$build_root/render-vignette.R"
status_before="$build_root/status-before"
status_after="$build_root/status-after"

cleanup() {
    case "$build_root" in
        /tmp/probatch-build-vignettes.*)
            if [[ -d "$build_root" && ! -L "$build_root" ]]; then
                rm -rf -- "$build_root"
            fi
            ;;
        *)
            printf 'Refusing to clean unexpected build path: %s\n' \
                "$build_root" >&2
            ;;
    esac
}
trap cleanup EXIT

install -d -m 700 -- "$stage_root"
for package_input in DESCRIPTION NAMESPACE R data inst vignettes; do
    if [[ -e "$package_input" ]]; then
        cp -a -- "$package_input" "$stage_root/"
    fi
done

mapfile -t staged_vignettes < <(
    find "$stage_root/vignettes" -maxdepth 1 -type f -name '*.Rmd' -print |
        LC_ALL=C sort
)
staged_vignette_count="$(
    printf '%s\n' "${staged_vignettes[@]}" |
        wc -l |
        tr -d ' '
)"
if (( staged_vignette_count == 0 )); then
    printf 'No staged vignette sources were found.\n' >&2
    exit 1
fi

cat > "$render_script" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Expected stage root, vignette source, and output directory")
}
stage_root <- normalizePath(args[[1L]], mustWork = TRUE)
source <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
intermediates_dir <- file.path(output_dir, "intermediates")
dir.create(intermediates_dir, recursive = TRUE, showWarnings = FALSE)

pkgload::load_all(
    stage_root,
    compile = FALSE,
    attach = TRUE,
    export_all = FALSE,
    export_imports = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE,
    quiet = TRUE,
    warn_conflicts = TRUE,
    debug = FALSE
)
options(mc.cores = 1L)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

rendered <- rmarkdown::render(
    input = source,
    output_format = NULL,
    output_dir = output_dir,
    intermediates_dir = intermediates_dir,
    knit_root_dir = dirname(source),
    envir = new.env(parent = globalenv()),
    clean = TRUE,
    quiet = TRUE
)
rendered <- normalizePath(rendered, mustWork = TRUE)
output_root <- paste0(normalizePath(output_dir, mustWork = TRUE), .Platform$file.sep)
if (!startsWith(rendered, output_root)) {
    stop("Rendered output escaped the isolated output directory: ", rendered)
}
cat(sprintf(
    "BUILT\t%s\t%s\t%s\n",
    basename(source),
    basename(rendered),
    unname(tools::md5sum(rendered))
))
RSCRIPT

safe_status_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
    "${safe_status_paths[@]}" > "$status_before"
pixi_before="$(sha256sum pixi.toml pixi.lock)"

for staged_vignette in "${staged_vignettes[@]}"; do
    vignette_name="$(basename "$staged_vignette" .Rmd)"
    output_dir="$build_root/output/$vignette_name"
    pixi run --as-is Rscript --vanilla \
        "$render_script" \
        "$stage_root" \
        "$staged_vignette" \
        "$output_dir"
done

git --no-optional-locks status --porcelain=v2 -z --untracked-files=all -- \
    "${safe_status_paths[@]}" > "$status_after"
cmp -s -- "$status_before" "$status_after" || {
    printf 'Vignette rendering changed the non-generated worktree.\n' >&2
    exit 1
}
pixi_after="$(sha256sum pixi.toml pixi.lock)"
if [[ "$pixi_before" != "$pixi_after" ]]; then
    printf 'Vignette rendering changed Pixi metadata.\n' >&2
    exit 1
fi

printf 'Built %d vignette sources in separate clean R processes without changing the source tree.\n' \
    "$staged_vignette_count"
```

## Diagnose Build Failure

```toml donna
id = "diagnose_build_failure"
kind = "donna.lib.request_action"
```

The isolated vignette build failed.

Standard output:

```text
{{ donna.lib.task_variable("build_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("build_stderr") }}
```

1. Identify the failing vignette, chunk, API, dependency, reference, or supporting asset.
2. Query all configured relations separately for every artifact before editing it.
3. Repair only authorized vignette code, prose, assets, package sources, package data, dependency metadata, or this workflow. Preserve unrelated user changes and the Git index.
4. For an R implementation change, preserve Roxygen2 and registration conventions in `R/`, add or update focused testthat coverage, and use the required test and documentation-source workflows before returning here.
5. Do not inspect generated manuals or `NAMESPACE`, generate documentation, install packages, use the network, or run package-wide checks.
6. After a repair, `{{ donna.lib.goto("verify_vignette_environment") }}`.
7. If the failure requires a maintainer decision, generated documentation, unavailable dependency installation, or work outside the authorized scope, `{{ donna.lib.goto("record_blocker") }}`.

## Review Build Results

```toml donna
id = "review_build_results"
kind = "donna.lib.request_action"
```

The isolated build reported:

```text
{{ donna.lib.task_variable("build_stdout") }}
```

1. Confirm every declared `.Rmd` source produced an isolated rendered output and no reference or supporting-asset error was suppressed.
2. Inspect the exact non-generated in-scope diff and any non-ignored untracked source artifacts without modifying the Git index. Preserve unrelated user changes.
3. Confirm the build used only declared dependencies, editable package sources, local data and assets, and no network or package installation.
4. If the results are complete and the current source changes are intentional, `{{ donna.lib.goto("verify_runtime_state") }}`.
5. If an in-scope repair is needed, make it after querying Depmesh and `{{ donna.lib.goto("verify_vignette_environment") }}`.
6. If the result cannot be accepted without maintainer action, `{{ donna.lib.goto("record_blocker") }}`.

## Verify Runtime State

```toml donna
id = "verify_runtime_state"
kind = "donna.lib.run_script"
save_stdout_to = "runtime_stdout"
save_stderr_to = "runtime_stderr"
goto_on_success = "final_review"
goto_on_failure = "repair_runtime_state"
timeout = 300
```

```bash donna script
#!/usr/bin/env bash
set -euo pipefail

generated_paths="$(
    find vignettes -mindepth 1 \
        \( \
            -type f \
            \( -name '*.html' -o -name '*.R' -o -name '.Rhistory' \) \
            -o -type d \
            \( -name '*_cache' -o -name '*_files' \) \
        \) -print
)"
if [[ -n "$generated_paths" ]]; then
    printf 'Generated vignette artifacts remain in the source tree:\n%s\n' \
        "$generated_paths" >&2
    exit 1
fi

for graphics_artifact in Rplots.pdf tests/testthat/Rplots.pdf; do
    if [[ -e "$graphics_artifact" || -L "$graphics_artifact" ]]; then
        printf 'Unexpected graphics artifact remains: %s\n' \
            "$graphics_artifact" >&2
        exit 1
    fi
done

safe_diff_paths=(
    .
    ':(top,exclude)man/**'
    ':(top,exclude)NAMESPACE'
)
git --no-optional-locks diff --check -- "${safe_diff_paths[@]}"
git --no-optional-locks diff --cached --check -- "${safe_diff_paths[@]}"

nonignored_vignette_paths="$(
    git --no-optional-locks ls-files --others --exclude-standard -- vignettes
)"
if printf '%s\n' "$nonignored_vignette_paths" |
    rg '(^|/)([^/]+(_cache|_files)(/|$)|[^/]+[.](html|R)$|[.]Rhistory$)'; then
    printf 'Non-ignored generated vignette artifacts remain.\n' >&2
    exit 1
fi

printf 'No generated vignette, graphics, cache, or HTML artifacts remain; non-generated diff hygiene passes.\n'
```

## Repair Runtime State

```toml donna
id = "repair_runtime_state"
kind = "donna.lib.request_action"
```

Final runtime-state verification failed.

Standard output:

```text
{{ donna.lib.task_variable("runtime_stdout") }}
```

Standard error:

```text
{{ donna.lib.task_variable("runtime_stderr") }}
```

1. Identify whether the named path was created by this workflow or pre-existed it. Do not delete or overwrite an artifact unless its ownership is exact and removal is authorized.
2. Query Depmesh before editing any governed source. Preserve unrelated changes and the Git index.
3. Repair the source or workflow so future rendering remains isolated, then `{{ donna.lib.goto("verify_vignette_environment") }}`.
4. If safe cleanup or repair requires developer direction, `{{ donna.lib.goto("record_blocker") }}`.

## Final Review

```toml donna
id = "final_review"
kind = "donna.lib.request_action"
```

1. Review the complete representative result and exact in-scope source diff without modifying the Git index.
2. Confirm both vignette sources are executable and current with the public API, their references and assets resolved, and generated outputs remained isolated from the repository.
3. Confirm every repair was focused, governed, and verified by any additionally triggered project workflows.
4. If the workflow completion outcome is satisfied, `{{ donna.lib.goto("finish") }}`.
5. If a repairable discrepancy remains, `{{ donna.lib.goto("verify_vignette_environment") }}`.
6. If a maintainer decision is required, `{{ donna.lib.goto("record_blocker") }}`.

## Finish

```toml donna
id = "finish"
kind = "donna.lib.finish"
```

Every declared package vignette built successfully in an isolated clean R process using declared local dependencies, references and assets resolved, and no generated rendering artifact entered the source tree.

If this was the first representative execution and the catalog still says `implementation in progress`, mark it `implemented` only after this finish, then rerun workflow validation, view rendering, artifact listing, bidirectional governance queries, triggered specification and file-relation verification, and scoped diff checks before reporting completion.

## Record Blocker

```toml donna
id = "record_blocker"
kind = "donna.lib.request_action"
```

1. Record any unresolved repository finding that no remaining workflow resolves in `{{ donna.lib.path("@/workflows/reports/vignette-open-questions.md") }}`. Query Depmesh before creating or editing that report, keep it non-executable, and do not add it to the workflow catalog.
2. Do not record transient command output as a repository question. Preserve exact evidence in the final report instead.
3. Then `{{ donna.lib.goto("blocked") }}`.

## Blocked

```toml donna
id = "blocked"
kind = "donna.lib.finish"
```

The vignette-build workflow cannot continue without maintainer input. Report completed gates, the exact dependency, source, asset, specification, or environment blocker, affected editable artifacts, and the decision required. Do not inspect generated manuals, modify `NAMESPACE`, install packages, use the network, stage files, create commits, or reset the Donna session.
