# Data dictionary

## Analysis units

Each record represents one pathogen-by-treatment observation. `Site` and
`Pathogen` jointly identify records when model predictions are joined to the
updated genomic-LRV labels.

## Core variables

| Variable | Definition |
|---|---|
| `Site` | Sampling site or treatment-pair identifier |
| `Pathogen` | Pathogen taxon |
| `Disinfection_Type` | UV, Cl, O3, UV_Cl, O3_Cl, or UF |
| `breadth_inf` | Initial genome breadth in the influent |
| `Delta_Breadth` | Change in genome breadth across treatment |
| `Delta_SNV_Density` | Change in SNV density across treatment |
| `Delta_Pi` | Change in nucleotide diversity across treatment |
| `LRV_point_new` | Point-estimate genomic log-removal value used as the model outcome |
| `LRV_status` | Label status; only `Point estimate` records are used for model fitting and scored validation |

## Input files

### `data/processed/xgboost_label_audit_all_447.csv`

Audited genomic-LRV labels and predictor fields for all 447 candidate training
records. The training scripts retain 430 point-estimate labels.

### `data/processed/training_full_model_completecase_243.csv`

The 243 point-estimate records with complete values for all full-model
predictors.

### `data/processed/training_reduced_model_430.csv`

All 430 point-estimate records used by the reduced model after excluding
`Delta_Pi`.

### `data/raw/data_no_internal_standard.xlsx`

Blind-test predictor data. The strict full-model analysis evaluates 157 complete
records. The imputation-assisted full model and the reduced model each evaluate
300 records with point-estimate outcomes.

## Output groups

### Training and internal validation

- `training_fit_predictions.csv`: apparent fitted values.
- `ten_fold_cv_predictions.csv`: pooled out-of-fold predictions.
- `training_and_cv_metrics.csv`: Pearson r², RMSE, MAE, and bias.
- `feature_importance_raw.csv`: encoded-feature XGBoost gain.
- `feature_importance_grouped.csv`: gain aggregated to manuscript predictors.

### Blind validation

- `blind_primary_no_imputation_predictions.csv`: primary strict validation.
- `blind_sensitivity_imputed_predictions.csv`: imputation-assisted validation.
- `blind_reduced_model_predictions.csv`: reduced-model validation.
- `blind_metrics_available_samples.csv`: overall and treatment-specific metrics.
- `blind_metrics_common_complete_subset.csv`: direct model comparison on the
  same 157 complete records.

### Missing-data audit

- `blind_imputation_audit_long.csv`: record- and feature-level imputation source.
- `blind_imputation_summary_all_features.csv`: counts by imputation source.
- `blind_imputation_summary_evaluated.csv`: complete versus imputed coverage.
- `imputation_assisted_metrics_by_input_status.csv`: performance reported
  separately for complete and imputed records.

### Model interpretation and error analysis

- `SHAP_full_model_243_source_data.csv`: row-level grouped TreeSHAP values.
- `SHAP_full_model_243_importance.csv`: mean absolute SHAP summaries.
- `Prediction_Error_Strict_No_Imputation_Source_Data.csv`: row-level absolute
  errors for the strict blind test.
- `Prediction_Error_Strict_No_Imputation_Group_Summary.csv`: group sample sizes,
  mean absolute errors, and medians.
- `Prediction_Error_Strict_No_Imputation_Pairwise_Wilcoxon.csv`: pairwise
  Wilcoxon tests with Benjamini–Hochberg correction.

## Interpretation boundary

`Observed_LRV` and `Predicted_LRV` refer to DNA-based genomic removal. They do
not directly represent viable-cell inactivation or infectivity reduction.

In the core CSV metric tables, `R2` is squared Pearson correlation between
observed and predicted values. The optional SI export additionally reports
predictive R², calculated as \(1-\mathrm{SSE}/\mathrm{SST}\).
