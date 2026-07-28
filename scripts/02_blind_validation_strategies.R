rm(list = ls())
source("scripts/model_utils.R")

training_full <- read_csv(
  "data/processed/training_full_model_completecase_243.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
blind_original <- read_excel(
  "data/raw/data_no_internal_standard.xlsx"
) %>%
  normalize_categories() %>%
  mutate(.blind_row_id = row_number())
labels <- read_csv(
  "data/processed/xgboost_label_audit_all_447.csv",
  show_col_types = FALSE
) %>%
  filter(LRV_status == "Point estimate") %>%
  transmute(
    Site,
    Pathogen,
    Observed_LRV = LRV_point_new,
    LRV_status
  )

model_full <- xgb.load(
  "models/XGBoost_full_completecase_243.ubj"
)
model_reduced <- xgb.load(
  "models/XGBoost_reduced_no_DeltaPi_430.ubj"
)

# Strategy 1: primary strict external validation with no imputation.
blind_primary_features <- blind_original %>%
  filter(
    complete.cases(
      across(all_of(required_columns(FULL_FEATURES)))
    )
  )
blind_primary_features$Predicted_LRV <- predict_xgboost(
  model_full,
  blind_primary_features,
  FULL_FEATURES
)
blind_primary <- blind_primary_features %>%
  inner_join(labels, by = c("Site", "Pathogen")) %>%
  filter(is.finite(Observed_LRV), is.finite(Predicted_LRV)) %>%
  mutate(
    Strategy = "Strict no-imputation",
    Input_status = "Complete"
  )

# Strategy 2: imputation-assisted sensitivity validation.
imputed <- impute_blind_from_training(
  blind_original %>% select(-.blind_row_id),
  training_full,
  FULL_FEATURES
)
blind_imputed_features <- imputed$data %>%
  mutate(
    .blind_row_id = row_number(),
    Input_status = if_else(
      Delta_Pi_missing_original == 1,
      "Imputed",
      "Complete"
    )
  )
blind_imputed_features$Predicted_LRV <- predict_xgboost(
  model_full,
  blind_imputed_features,
  FULL_FEATURES
)
blind_imputed <- blind_imputed_features %>%
  inner_join(labels, by = c("Site", "Pathogen")) %>%
  filter(is.finite(Observed_LRV), is.finite(Predicted_LRV)) %>%
  mutate(Strategy = "Imputation-assisted")

# Strategy 3: deployment-oriented model trained on all 430 point labels and
# excluding Delta_Pi. No blind-test imputation is required.
assert_complete_features(
  blind_original,
  REDUCED_FEATURES,
  "reduced-model blind data"
)
blind_reduced_features <- blind_original
blind_reduced_features$Predicted_LRV <- predict_xgboost(
  model_reduced,
  blind_reduced_features,
  REDUCED_FEATURES
)
blind_reduced <- blind_reduced_features %>%
  inner_join(labels, by = c("Site", "Pathogen")) %>%
  filter(is.finite(Observed_LRV), is.finite(Predicted_LRV)) %>%
  mutate(
    Strategy = "Reduced model (no Delta_Pi)",
    Input_status = "No imputation required"
  )

metrics_available <- bind_rows(
  metrics_by_disinfection(
    blind_primary,
    "Full model (243)",
    "Strict no-imputation"
  ),
  metrics_by_disinfection(
    blind_imputed,
    "Full model (243)",
    "Imputation-assisted"
  ),
  metrics_by_disinfection(
    blind_reduced,
    "Reduced model (430)",
    "Reduced model (no Delta_Pi)"
  )
)

# Direct model comparison on the same complete-case records.
common_ids <- blind_primary %>%
  select(.blind_row_id)
blind_imputed_common <- blind_imputed %>%
  semi_join(common_ids, by = ".blind_row_id")
blind_reduced_common <- blind_reduced %>%
  semi_join(common_ids, by = ".blind_row_id")
metrics_common <- bind_rows(
  metrics_by_disinfection(
    blind_primary,
    "Full model (243)",
    "Strict no-imputation",
    "Common complete-case subset"
  ),
  metrics_by_disinfection(
    blind_imputed_common,
    "Full model (243)",
    "Imputation-assisted",
    "Common complete-case subset"
  ),
  metrics_by_disinfection(
    blind_reduced_common,
    "Reduced model (430)",
    "Reduced model (no Delta_Pi)",
    "Common complete-case subset"
  )
)

sensitivity_status_metrics <- bind_rows(lapply(
  c("Complete", "Imputed"),
  function(status) {
    d <- blind_imputed %>% filter(Input_status == status)
    metrics_by_disinfection(
      d,
      "Full model (243)",
      paste0("Imputation-assisted: ", status)
    )
  }
))

imputation_audit <- imputed$audit
imputation_summary_all <- imputation_audit %>%
  filter(Was_missing) %>%
  count(Feature, Imputation_source, name = "n_imputed") %>%
  mutate(
    Denominator = nrow(blind_original),
    Imputed_pct = n_imputed / Denominator * 100
  )
imputation_summary_evaluated <- blind_imputed %>%
  count(Input_status, name = "n_evaluated") %>%
  mutate(
    Denominator = nrow(blind_imputed),
    Percentage = n_evaluated / Denominator * 100
  )

blind_all <- bind_rows(
  blind_primary,
  blind_imputed,
  blind_reduced
)
all_metrics <- bind_rows(metrics_available, metrics_common)

write_csv(
  blind_primary,
  "results/tables/blind_primary_no_imputation_predictions.csv"
)
write_csv(
  blind_imputed,
  "results/tables/blind_sensitivity_imputed_predictions.csv"
)
write_csv(
  blind_reduced,
  "results/tables/blind_reduced_model_predictions.csv"
)
write_csv(
  blind_all,
  "results/tables/blind_all_strategies_predictions.csv"
)
write_csv(
  metrics_available,
  "results/tables/blind_metrics_available_samples.csv"
)
write_csv(
  metrics_common,
  "results/tables/blind_metrics_common_complete_subset.csv"
)
write_csv(
  all_metrics,
  "results/tables/blind_metrics_all_comparisons.csv"
)
write_csv(
  sensitivity_status_metrics,
  "results/tables/imputation_assisted_metrics_by_input_status.csv"
)
write_csv(
  imputation_audit,
  "results/tables/blind_imputation_audit_long.csv"
)
write_csv(
  imputation_summary_all,
  "results/tables/blind_imputation_summary_all_features.csv"
)
write_csv(
  imputation_summary_evaluated,
  "results/tables/blind_imputation_summary_evaluated.csv"
)

print(imputation_summary_all)
print(imputation_summary_evaluated)
print(metrics_available %>% filter(Group == "Overall"))
print(sensitivity_status_metrics %>% filter(Group == "Overall"))
