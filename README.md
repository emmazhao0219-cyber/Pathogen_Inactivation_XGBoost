# Cross-species XGBoost prediction of DNA-based genomic removal

This repository contains the data, code, trained models, validation outputs, and
figure source data for a cross-species XGBoost analysis of pathogen-associated
**genomic log-removal values (genomic LRVs)** in wastewater treatment.

The models use metagenome-derived damage signatures, particularly changes in
genome breadth, SNV density, and nucleotide diversity, together with initial
genome breadth and disinfection type. The response is DNA-based genomic removal.
It should not be interpreted as viable-cell inactivation, infectivity loss, or a
direct health-risk endpoint.

## Analysis design

| Analysis | Training data | Blind-test coverage | Missing-data handling | Intended role |
|---|---:|---:|---|---|
| Full model, strict validation | 243 complete records | 157 complete records | None | Primary external validation |
| Full model, imputation-assisted validation | Same 243 records | 300 records | Training-derived hierarchical medians | Sensitivity analysis |
| Reduced model without `Delta_Pi` | 430 point-estimate records | 300 records | None | Coverage and applicability analysis |

The hierarchical imputation order is:

1. disinfection type × pathogen training-set median;
2. pathogen training-set median;
3. overall training-set median.

Missingness indicators are retained for audit and stratified reporting, but are
not supplied to the XGBoost model.

## Repository structure

```text
.
├── data/
│   ├── raw/                  # blind-test feature input
│   └── processed/            # analysis labels and model-ready datasets
├── models/                   # fitted full and reduced XGBoost models
├── results/
│   ├── figures/              # publication figures
│   ├── tables/               # predictions, metrics, SHAP and audit tables
│   └── si_internal_validation/  # optional SI export
├── scripts/
│   ├── 00_check_environment.R
│   ├── 01_model_training.R
│   ├── 02_blind_validation_strategies.R
│   ├── 03_data_visualization.R
│   ├── 04_shap_full_model.R
│   ├── 05_prediction_error_by_disinfection.R
│   ├── 06_si_internal_validation_material.R
│   └── model_utils.R
├── DATA_DICTIONARY.md
└── run_all.R
```

## Requirements

The workflow was tested with R 4.4.3. Required packages are declared in
[`DESCRIPTION`](DESCRIPTION):

```r
install.packages(c(
  "readr", "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "xgboost", "ggplot2", "svglite", "ragg", "ggbeeswarm",
  "systemfonts"
))
```

Arial is used when available. On systems without Arial, the scripts fall back to
the generic sans-serif family without changing the analysis.

## Reproduce the analysis

Clone the repository and run all commands from any location:

```bash
git clone https://github.com/emmazhao0219-cyber/Pathogen_Inactivation_XGBoost.git
cd Pathogen_Inactivation_XGBoost
Rscript scripts/00_check_environment.R
Rscript run_all.R
```

The default pipeline trains both models, performs stratified 10-fold
cross-validation, runs all three blind-test strategies, and regenerates the PNG,
PDF, and SVG outputs.

To additionally export 600-dpi TIFF figures:

```bash
EXPORT_TIFF=1 Rscript run_all.R
```

To also create the SI internal-validation package:

```bash
BUILD_SI=1 Rscript run_all.R
```

Both options can be enabled together.

## Reproducibility settings

- Random seed: `2026`.
- XGBoost objective: `reg:squarederror`.
- Hyperparameters: `eta = 0.1`, `max_depth = 4`,
  `min_child_weight = 3`, `subsample = 0.8`.
- Boosting rounds: `100`.
- Internal validation: disinfection-stratified 10-fold cross-validation.
- Disinfection type: one-hot encoded with UV as the reference level.
- Model outcome: point-estimate genomic LRV only.

## Key validation results

| Model and evaluation | n | Pearson r² | RMSE | MAE | Bias |
|---|---:|---:|---:|---:|---:|
| Full model, apparent training fit | 243 | 0.942 | 0.125 | 0.091 | −0.001 |
| Full model, stratified 10-fold CV | 243 | 0.454 | 0.364 | 0.242 | 0.003 |
| Full model, strict blind validation | 157 | 0.511 | 0.268 | 0.207 | −0.046 |
| Full model, imputation-assisted blind validation | 300 | 0.339 | 0.315 | 0.232 | −0.080 |
| Reduced model, blind validation | 300 | 0.335 | 0.307 | 0.225 | −0.011 |

The apparent training fit is not a measure of generalization. The strict
no-imputation blind test is the primary external validation. Imputation-assisted
results are reported separately because the 143 imputed records were less
predictable than the 157 complete records.

In the core result tables, `R2` denotes squared Pearson correlation. Predictive
R², together with Pearson r², is reported in the optional SI internal-validation
export.

See [`RESULTS_SUMMARY.md`](RESULTS_SUMMARY.md) for interpretation and
[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) for input and output definitions.

## Main outputs

- `results/tables/training_and_cv_metrics.csv`
- `results/tables/blind_metrics_available_samples.csv`
- `results/tables/blind_metrics_common_complete_subset.csv`
- `results/tables/imputation_assisted_metrics_by_input_status.csv`
- `results/tables/blind_imputation_audit_long.csv`
- `results/tables/SHAP_full_model_243_source_data.csv`
- `results/figures/Blind_Primary_No_Imputation_Faceted.*`
- `results/figures/SHAP_Summary_Full_Model_243.*`
- `results/figures/Prediction_Error_By_Disinfection_Strict_No_Imputation.*`

## Citation

If you use this repository, please cite the associated manuscript and the
software metadata in [`CITATION.cff`](CITATION.cff). The manuscript DOI can be
added to `CITATION.cff` after publication.
