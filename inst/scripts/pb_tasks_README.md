# pb_tasks.yaml — Combination Reference

This file documents all batch-correction combinations defined in `pb_tasks.yaml`.

---

## Workflow Usage

In `inst/scripts/batch_correction_workflow.R`, the default configuration points to this file:

```r
correction_tasks_yaml <- "inst/scripts/pb_tasks.yaml"
correction_task_labels <- NULL
```

With `correction_task_labels <- NULL`, all expanded tasks from this YAML are used.

---

## Structure Overview

The task grid expands into three branches:

1. **Direct correction** — correction methods run directly on `log2_on_raw`
2. **Imputation-only** — `log2_on_raw` is imputed (no correction)
3. **Imputation + correction** — `log2_on_raw` is first imputed, then corrected

---

## Profiles (preprocessing anchors)

| ID | Step | Behavior |
|----|------|----------|
| `log2_on_raw` | *(none)* | Pass-through anchor; uses `input_assay` directly |

---

## Imputation Methods

| ID | Step | Notes |
|----|------|-------|
| `omicsGMFImpute` | omicsGMFImpute | GMF-based; uses `batch_col` design formula |
| `PRONEImpute_conditioned` | PRONEImpute | Conditioned on `condition_col` |
| `PRONEImpute_global` | PRONEImpute | Global mode (`condition_col = null`) |
| `MFimpute` | MFimpute | Random forest imputation via ranger |

---

## Correction Methods

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

---

## Branch 1 — Direct Correction

**1 profile × 12 methods = 12 max combinations**

| Profile | Correction Methods |
|---------|--------------------|
| `log2_on_raw` | all 12 above |

---

## Branch 2 — Imputation Only

**1 profile × 4 imputation methods = 4 combinations**
(only when `combinations.include_imputation_outputs: true`)

| Profile | Imputation Method |
|---------|-------------------|
| `log2_on_raw` | omicsGMFImpute |
| `log2_on_raw` | PRONEImpute_conditioned |
| `log2_on_raw` | PRONEImpute_global |
| `log2_on_raw` | MFimpute |

---

## Branch 3 — Imputation + Correction

**1 profile × 4 imputation methods × 12 correction methods = 48 max combinations**

| Profile | Imputation | Correction Methods |
|---------|------------|--------------------|
| `log2_on_raw` | omicsGMFImpute | all 12 above |
| `log2_on_raw` | PRONEImpute_conditioned | all 12 above |
| `log2_on_raw` | PRONEImpute_global | all 12 above |
| `log2_on_raw` | MFimpute | all 12 above |

---

## Grand Total

| Branch | Max | Runnable (with defaults) |
|--------|-----|--------------------------|
| Direct correction | 12 | 8 |
| Imputation-only | 4 | 4 |
| Imputation + Correction | 48 | 32 |
| **Total** | **64** | **44** |

Runnable with defaults excludes methods requiring unset params (`mComBat`, `BERT_ref`, `RUVIIIC_k5`, `RUVIIIC_k3`).

---

## Enabling Skipped Methods

| Method | Parameter to set in `user_params` |
|--------|-----------------------------------|
| `mComBat` | `mComBat_center` |
| `BERT_ref` | `bert_reference_name` |
| `RUVIIIC_k5`, `RUVIIIC_k3` | `replicate_col` + `negative_control_features` |
