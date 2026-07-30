# Standalone Core baseline

## Goal of the document

This document establishes the observable package, ownership, data, extension, embedding, and dependency contracts for a standalone `proBatch` Core that can be consumed in one direction by `proBatchBench`.

## Scope

The scope of this specification includes the Core-to-Bench ownership boundary; matrix, sample-alignment, assay-name, lineage, and missing-value invariants; the provider-neutral extension surface required by the pinned Bench consumer; Core-owned PCA, t-SNE, and UMAP diagnostics; and the editable package artifacts and dependencies that expose those contracts.

Provider implementations, benchmark orchestration and metrics, release version selection, generated files under `./man/`, generated `./NAMESPACE`, and incidental private-helper structure are out of scope.

## Dictionary

- `Core` — the independently installable `proBatch` package in this repository.
- `Bench` — the downstream `proBatchBench` package, which may consume public Core interfaces.
- `provider` — a package that registers a transformation implemented outside Core.
- `stored assay` — an assay materialized in a `ProBatchFeatures` object.
- `virtual target` — an operation-log result that can be replayed but is not currently materialized as an assay.
- `stable lineage origin` — the operation fields `from`, `to`, `step`, `fun`, `pkg`, and `params`; execution timestamps are not stable lineage origin fields.
- `structured step result` — transformed data paired with provider-neutral named or unnamed artifacts that are not part of the assay matrix.

## Package and ownership boundary

`proBatch` MUST install, load, and provide its Core behavior without `proBatchBench` or a provider-specific engine.

Bench MAY import public Core interfaces. Core MUST NOT import Bench, call Bench code, or declare Bench or provider-specific correction, imputation, normalization, or benchmark engines as package dependencies.

Core MUST own general containers, representation conversion, transformations, correction primitives, missing-value policies, lineage, sample alignment, diagnostic plots, and PCA, t-SNE, and UMAP over matrices and `ProBatchFeatures`.

Bench MUST own provider integrations, provider-specific dependencies, benchmark planning and execution, evaluation metrics, and downstream presentation that is not a general Core diagnostic.

An extension point added for Bench MUST be provider-neutral, usable with an in-package fake provider, and justified by an observed downstream call. Core MUST NOT copy provider implementations or orchestration from Bench.

## Data and lineage contracts

### Matrix orientation and sample alignment

A matrix accepted or returned by a Core transformation or embedding API MUST represent features in rows and samples in columns.

Feature and sample identifiers required by an API MUST be non-missing, non-empty, and unique. An API MUST NOT align sample data by position when identifiers are available or required.

Sample annotation MUST be aligned to matrix columns by `sample_id_col`. An API that explicitly documents row-name fallback MAY use unique, non-empty annotation row names when `sample_id_col` is absent. Duplicate annotation identifiers or a matrix sample without annotation MUST cause an error. Extra annotation rows MAY be ignored, and accepted annotation MUST be reordered to matrix sample order.

For `ProBatchFeatures` diagnostics, `pbf_name = NULL` MUST select assays in object order. An explicit assay selection MUST preserve first-occurrence request order. A single selected assay MUST return the single-assay result; a multi-assay result MUST preserve assay names and order.

### Assay names and target collisions

An explicit assay or result name MUST be one non-missing, non-empty character value.

Stored assay names and virtual targets share one result namespace. A new explicit name that is already stored or reserved by the operation log MUST cause an error unless the request is an exact idempotent retry.

An existing-target retry is idempotent only when the requested data, parent assay, and stable lineage origin are identical to the existing target. Conflicting data, a different parent, or different stable origin fields MUST cause an error.

An API MAY disambiguate a generated name when that behavior is documented. It MUST NOT silently suffix, rewrite, or otherwise disambiguate a caller-supplied explicit final name. This rule applies equally to `pb_filterNA()` and `pb_groupfilterNA()`.

### Operation lineage

Every non-self result in the operation log MUST have at most one non-self parent.

Lineage-edge identity MUST compare every stable lineage origin field. Recording an identical edge MAY be idempotent; recording a second origin for the same result MUST cause an error.

Self edges MAY record in-place operations without changing the non-self parent. Stored and virtual operations MUST obey the same conflict rules.

`get_chain()`, `pb_pipeline_name()`, and virtual replay through `pb_assay_matrix()` MUST follow one unambiguous lineage. Cyclic or conflicting lineage MUST cause an error rather than choosing an edge by row order.

### Missing values

Provider-neutral transformation and correction interfaces MUST use the canonical policies `error`, `keep`, `drop_features`, and `fill`.

`error` MUST reject missing input, `keep` MUST pass missing values to a method that supports them, `drop_features` MUST remove incomplete feature rows while preserving retained order, and `fill` MUST replace missing values with an explicit numeric `fill_value`. Removing every feature MUST cause an error.

PCA, t-SNE, and UMAP MUST share the plotting-family `fill_the_missing` behavior. Their default MUST remain `-1`; a numeric scalar fills missing values, while `NULL` or `FALSE` removes incomplete features. Each embedding MUST warn before applying non-error missing-value handling.

## Provider-neutral extension surface

### Step registration lifecycle

`pb_register_step()` MUST accept a canonical name and function and MAY record a provider package, step kind, required packages, aliases, and display label. Omitting provider metadata MUST remain compatible with existing two-argument registrations by inferring a provider without loading another package.

A registration MUST NOT install, load, or probe a provider's optional engines. Exact repeated registration MUST be idempotent. A canonical-name or alias collision with a different registration MUST cause an error.

`pb_unregister_steps(package)` MUST remove all canonical registrations and aliases owned by that provider and MUST be idempotent when none exist. Provider `.onLoad` and `.onUnload` hooks MUST therefore be able to register and clean up one provider without changing Core or another provider's registrations.

`pb_list_steps()` MUST return canonical names by default. With `details = TRUE`, it MUST return, at minimum, canonical `name`, provider `package`, `kind`, `label`, `requires`, `aliases`, and current `available` state. Pattern filtering MUST match canonical names, aliases, and labels. Availability filtering MUST NOT mutate the registry.

`pb_has_step()` MUST resolve canonical names and aliases and MUST optionally distinguish registration from current availability.

A provider registration is available only when its provider is loaded and every declared requirement is available; Core and ordinary user registrations MAY be available without a separate provider namespace. Invocation MUST report an unavailable provider or requirement before running the function.

Invoking an alias through `pb_transform()` or `pb_eval()` MUST record the canonical function name and provider package. Replay MUST resolve that recorded provider and MUST fail with actionable guidance when the provider is no longer registered or available.

### Matrix-method adaptation

`pb_apply_matrix_method()` MUST accept a numeric feature-by-sample matrix or a long data frame. A `ProBatchFeatures` object MUST instead use `pb_transform()` so storage and lineage are preserved.

Long input MUST identify feature, sample, and numeric measure columns. Duplicate feature/sample keys MUST cause an error rather than implicit aggregation.

The adapter MUST construct a numeric feature-by-sample matrix, apply the canonical missing-value policy, align sample annotation by identifier, call the method with the matrix and aligned annotation, and validate the returned data.

Method output MUST contain only known feature and sample identifiers and MUST preserve input-relative order; an ordered subset is allowed. A non-`keep` result MUST NOT introduce missing values.

For matrix input, the adapter MUST return a matrix. For long input, it MUST restore rows without a join or Cartesian expansion, preserve retained rows' relative order, and support documented minimal, original-column, and non-conflicting annotation-column retention modes.

### Structured step results

`pb_step_result(data, artifacts = list())` MUST construct a validated structured step result, and `artifacts` MUST be a list.

Matrix adaptation MUST preserve a structured result while applying the same validation and representation restoration to its `data`.

When `pb_transform()` materializes a structured result, it MUST store the transformed matrix as the assay data and preserve the artifacts in that assay's metadata under `pb_step_artifacts`. Artifact persistence is not required for `pb_eval()` or an unmaterialized virtual result.

## Core embedding diagnostics

Core MUST export the `plot_TSNE()` and `plot_UMAP()` S3 generics and register `default` and `ProBatchFeatures` methods in generated namespace output after maintainer regeneration.

The default methods MUST treat matrix columns as observations, align annotation by sample identifiers, and preserve sample order in returned plot data.

`plot_TSNE()` MUST use `Rtsne` and `plot_UMAP()` MUST use `umap`. The backend packages MUST be optional package dependencies guarded with actionable errors. `plotly` MUST be required only when interactive rendering is requested.

Static rendering MUST return a `ggplot` by default. An explicit interactive-rendering option MUST return a `plotly` object for one assay. Multiple interactive assays MUST return a named ordered list unless subplot combination is explicitly requested.

For multiple static assays, `return_gridExtra = TRUE` MUST expose both the arranged grob and a named ordered plot list. For multiple interactive assays, `return_subplots = TRUE` MUST combine plots with documented column and shared-axis controls.

t-SNE MUST require at least two samples and at least two requested embedding dimensions. A perplexity above `max(1, floor((n_samples - 1) / 3))` MUST be reduced to that value with a warning. Backend arguments not owned by Core MUST be forwarded without changing matrix orientation.

UMAP MUST require at least two requested components and MUST forward documented neighborhood, distance, metric, spread, learning-rate, and backend arguments without changing matrix orientation.

t-SNE MUST expose an explicit `random_seed`, and UMAP MUST expose an explicit `random_state`. Repeated calls with the same data, arguments, seed, and backend versions MUST produce the same embedding coordinates. A missing seed MUST leave backend randomness under the caller's control.

Both embeddings MUST support factor-like and numeric color annotations, optional factor-like shapes, requested assay ordering, per-assay titles and arguments, and the plotting-family static and interactive color semantics.

## Dependency and editable-source contract

`./DESCRIPTION` MUST remain the authority for package name, version, supported R version, and R/Bioconductor dependencies. The rejected split's package version and R-version choices MUST NOT change those values.

`Rtsne`, `umap`, and `plotly` MUST be declared in `Suggests` and guarded at runtime. Direct Pixi dependencies MUST make their tested paths available without making them required to install Core.

The declared `testthat` minimum and the direct Pixi constraint MUST be at least `3.2.0`.

Core source and tests MUST NOT directly use `plyr` or `reshape2`, and those packages MUST NOT remain direct Core or Pixi requirements solely for removed uses. Dependency removal beyond these proven cases MUST require source-level evidence.

`./.Rbuildignore` MUST exclude `./.pixi`, `./pixi.toml`, and `./pixi.lock` from source-package archives.

Each public R symbol MUST have one effective top-level definition so source order does not select behavior. Existing exported correction compatibility wrappers MUST remain documented forwarding wrappers until a separately specified removal.

Public exports, methods, parameters, and examples MUST be maintained in Roxygen2 sources under `./R/`; agents MUST leave `./NAMESPACE` and `./man/` generation to the maintainer.

Changes to Core behavior MUST include focused testthat coverage and synchronized editable Roxygen2 sources, vignettes, dependency metadata, Depmesh relations, and release notes when those artifact families are affected.
