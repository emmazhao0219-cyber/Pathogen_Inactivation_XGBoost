rm(list = ls())
source("scripts/model_utils.R")

audit <- read_csv(
  "data/processed/xgboost_label_audit_all_447.csv",
  show_col_types = FALSE
)
point_estimates <- audit %>%
  filter(LRV_status == "Point estimate") %>%
  mutate(LRV_model = LRV_point_new) %>%
  normalize_categories()

training_full <- point_estimates %>%
  filter(
    complete.cases(
      across(all_of(required_columns(FULL_FEATURES)))
    )
  )
training_reduced <- point_estimates %>%
  filter(
    complete.cases(
      across(all_of(required_columns(REDUCED_FEATURES)))
    )
  )

stopifnot(
  nrow(point_estimates) == 430,
  nrow(training_full) == 243,
  nrow(training_reduced) == 430
)
assert_complete_features(
  training_full,
  FULL_FEATURES,
  "full-model training data"
)
assert_complete_features(
  training_reduced,
  REDUCED_FEATURES,
  "reduced-model training data"
)

write_csv(
  training_full,
  "data/processed/training_full_model_completecase_243.csv"
)
write_csv(
  training_reduced,
  "data/processed/training_reduced_model_430.csv"
)

model_full <- fit_xgboost(training_full, FULL_FEATURES)
model_reduced <- fit_xgboost(training_reduced, REDUCED_FEATURES)
xgb.save(
  model_full,
  "models/XGBoost_full_completecase_243.ubj"
)
xgb.save(
  model_reduced,
  "models/XGBoost_reduced_no_DeltaPi_430.ubj"
)

training_full$Predicted_LRV <- predict_xgboost(
  model_full,
  training_full,
  FULL_FEATURES
)
training_reduced$Predicted_LRV <- predict_xgboost(
  model_reduced,
  training_reduced,
  REDUCED_FEATURES
)

cv_full <- cross_validate_xgboost(
  training_full,
  FULL_FEATURES,
  "Full model (243)"
)
cv_reduced <- cross_validate_xgboost(
  training_reduced,
  REDUCED_FEATURES,
  "Reduced model (430)"
)

training_predictions <- bind_rows(
  training_full %>%
    transmute(
      Model = "Full model (243)",
      Site,
      Pathogen,
      Disinfection_Type,
      breadth_inf,
      Observed_LRV = LRV_model,
      Predicted_LRV
    ),
  training_reduced %>%
    transmute(
      Model = "Reduced model (430)",
      Site,
      Pathogen,
      Disinfection_Type,
      breadth_inf,
      Observed_LRV = LRV_model,
      Predicted_LRV
    )
)
cv_predictions <- bind_rows(cv_full, cv_reduced)

training_cv_metrics <- bind_rows(
  metric_row(
    training_full$LRV_model,
    training_full$Predicted_LRV,
    "Full model (243)",
    "Apparent training fit"
  ),
  metric_row(
    cv_full$Observed_LRV,
    cv_full$Predicted_LRV,
    "Full model (243)",
    "Stratified 10-fold CV"
  ),
  metric_row(
    training_reduced$LRV_model,
    training_reduced$Predicted_LRV,
    "Reduced model (430)",
    "Apparent training fit"
  ),
  metric_row(
    cv_reduced$Observed_LRV,
    cv_reduced$Predicted_LRV,
    "Reduced model (430)",
    "Stratified 10-fold CV"
  )
)

importance_full <- xgb.importance(
  feature_names = colnames(
    make_feature_matrix(training_full, FULL_FEATURES)
  ),
  model = model_full
) %>%
  as_tibble() %>%
  mutate(Model = "Full model (243)")
importance_reduced <- xgb.importance(
  feature_names = colnames(
    make_feature_matrix(training_reduced, REDUCED_FEATURES)
  ),
  model = model_reduced
) %>%
  as_tibble() %>%
  mutate(Model = "Reduced model (430)")
importance <- bind_rows(importance_full, importance_reduced) %>%
  group_by(Model) %>%
  mutate(Gain_pct = Gain / sum(Gain) * 100) %>%
  ungroup()

write_csv(
  training_predictions,
  "results/tables/training_fit_predictions.csv"
)
write_csv(
  cv_predictions,
  "results/tables/ten_fold_cv_predictions.csv"
)
write_csv(
  training_cv_metrics,
  "results/tables/training_and_cv_metrics.csv"
)
write_csv(
  importance,
  "results/tables/feature_importance_raw.csv"
)
print(training_cv_metrics)
