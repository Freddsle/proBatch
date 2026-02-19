# Batch Correction Workflow README

This README is the single reference for running/customizing the batch-correction workflow and for the `pb_tasks.yaml` task-grid specification in `inst/scripts`.

## File Roles

| File | Role | Typical edit frequency |
|------|------|------------------------|
| `inst/scripts/batch_correction_workflow_params.R` | User-facing run configuration (`workflow_params`) and entrypoint | Often (per run) |
| `inst/scripts/batch_correction_workflow.R` | Core workflow engine (diagnostics, correction, reporting) | Rarely |
| `inst/scripts/batch_correction_task_helpers.R` | Shared helpers for task expansion/execution and workflow utility functions | Rarely |
| `inst/scripts/pb_tasks.yaml` | Canonical task-grid template used by default | Rarely |

## Recommended Edit Policy

1. Edit run-specific settings in `inst/scripts/batch_correction_workflow_params.R`.
2. Keep `inst/scripts/batch_correction_workflow.R` and `inst/scripts/batch_correction_task_helpers.R` stable unless changing workflow logic.
3. Treat `inst/scripts/pb_tasks.yaml` as canonical. For run-specific task-graph changes, copy it and point `workflow_params$correction_tasks_yaml` to the copy.

## Quick Start

1. Set values in `workflow_params` inside `inst/scripts/batch_correction_workflow_params.R`.
2. Run:

```bash
Rscript inst/scripts/batch_correction_workflow_params.R
```

The params file resolves and sources `inst/scripts/batch_correction_workflow.R`, which applies overrides from `workflow_params`.

## Custom Task YAML Per Run

If you need run-specific task graph changes (profiles/combinations/method sets):

1. Copy:

```bash
cp inst/scripts/pb_tasks.yaml inst/scripts/pb_tasks.user.yaml
```

2. Edit `inst/scripts/pb_tasks.user.yaml`.
3. Set in `workflow_params`:

```r
correction_tasks_yaml = "inst/scripts/pb_tasks.user.yaml"
```

Use `correction_task_labels` only for exact label filtering of expanded tasks.

## What to Customize Where

- Use `workflow_params` for dataset paths, column names, plotting options, output folder, and task-file path.
- Use task YAML for method graph design (profiles, imputation/correction combinations, method-specific required params).
- If you only need to enable methods skipped by required params, set those required values in YAML `task_grid.user_params` (or in objects referenced by YAML).

## Task Grid Reference (`pb_tasks.yaml`)

### Workflow Usage

In `inst/scripts/batch_correction_workflow.R`, default task settings point to:

```r
correction_tasks_yaml <- "inst/scripts/pb_tasks.yaml"
correction_task_labels <- NULL
```

With `correction_task_labels <- NULL`, all expanded tasks from this YAML are used.

### Editability Policy

- `inst/scripts/pb_tasks.yaml` is the canonical template for this repository.
- For run-specific changes, copy this file (for example `inst/scripts/pb_tasks.user.yaml`) and edit the copy.
- Point `workflow_params$correction_tasks_yaml` in `inst/scripts/batch_correction_workflow_params.R` to your copy.
- Keep routine run settings (paths, columns, plotting, outputs) in `workflow_params`, not in this YAML.

### Structure Overview

The task grid expands into three branches:

1. **Direct correction** — correction methods run directly on `log2_on_raw`
2. **Imputation-only** — `log2_on_raw` is imputed (no correction)
3. **Imputation + correction** — `log2_on_raw` is first imputed, then corrected

### Profiles (preprocessing anchors)

| ID | Step | Behavior |
|----|------|----------|
| `log2_on_raw` | *(none)* | Pass-through anchor; uses `input_assay` directly |

### Imputation Methods

| ID | Step | Notes |
|----|------|-------|
| `omicsGMFImpute` | omicsGMFImpute | GMF-based; uses `batch_col` design formula |
| `PRONEImpute_conditioned` | PRONEImpute | Conditioned on `condition_col` |
| `PRONEImpute_global` | PRONEImpute | Global mode (`condition_col = null`) |
| `MFimpute` | MFimpute | Random forest imputation via ranger |

### Correction Methods

| ID | Step | Required params (skipped if null) |
|----|------|-----------------------------------|
| `ComBat` | combat | — |
| `mComBat` | mComBat | `mComBat_center` |
| `BERT_ComBat` | BERT | — |
| `BERT_limma` | BERT | — |
| `BERT_ref` | BERT | `bert_reference_name` |
| `limmaRBE` | limmaRBE | — |
| `RUVIIIC_k5` | RUVIIIC | `replicate_col`, `negative_control_features` |
| `RUVIIIC_k3` | RUVIIIC | `replicate_col`, `negative_control_features` |
| `NormAE` | NormAE | — |
| `PLSDAbatch` | PLSDAbatch | — |
| `sPLSDAbatch` | sPLSDAbatch | — |
| `omicsGMFcor` | omicsGMFcor | — |

Methods with missing required params are skipped when `settings.skip_invalid: true` (default).

### Branch 1 — Direct Correction

**1 profile × 12 methods = 12 max combinations**

| Profile | Correction Methods |
|---------|--------------------|
| `log2_on_raw` | all 12 above |

### Branch 2 — Imputation Only

**1 profile × 4 imputation methods = 4 combinations**
(only when `combinations.include_imputation_outputs: true`)

| Profile | Imputation Method |
|---------|-------------------|
| `log2_on_raw` | omicsGMFImpute |
| `log2_on_raw` | PRONEImpute_conditioned |
| `log2_on_raw` | PRONEImpute_global |
| `log2_on_raw` | MFimpute |

### Branch 3 — Imputation + Correction

**1 profile × 4 imputation methods × 12 correction methods = 48 max combinations**

| Profile | Imputation | Correction Methods |
|---------|------------|--------------------|
| `log2_on_raw` | omicsGMFImpute | all 12 above |
| `log2_on_raw` | PRONEImpute_conditioned | all 12 above |
| `log2_on_raw` | PRONEImpute_global | all 12 above |
| `log2_on_raw` | MFimpute | all 12 above |

### Grand Total

| Branch | Max | Runnable (with defaults) |
|--------|-----|--------------------------|
| Direct correction | 12 | 8 |
| Imputation-only | 4 | 4 |
| Imputation + Correction | 48 | 32 |
| **Total** | **64** | **44** |

Runnable with defaults excludes methods requiring unset params (`mComBat`, `BERT_ref`, `RUVIIIC_k5`, `RUVIIIC_k3`).

### Enabling Skipped Methods

| Method | Parameter to set in `user_params` |
|--------|-----------------------------------|
| `mComBat` | `mComBat_center` |
| `BERT_ref` | `bert_reference_name` |
| `RUVIIIC_k5`, `RUVIIIC_k3` | `replicate_col` + `negative_control_features` |

### Label Filtering Note

`correction_task_labels` uses exact task labels after expansion. If you provide labels that do not exist, no tasks are selected.
