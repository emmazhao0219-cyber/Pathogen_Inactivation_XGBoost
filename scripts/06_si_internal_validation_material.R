rm(list = ls())
source("scripts/model_utils.R")

output_dir <- file.path(
  "results",
  "si_internal_validation"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(output_dir, "figure"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(output_dir, "data"),
  recursive = TRUE,
  showWarnings = FALSE
)

training_source <- read_csv(
  "data/processed/training_full_model_completecase_243.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
training_predictions <- read_csv(
  "results/tables/training_fit_predictions.csv",
  show_col_types = FALSE
) %>%
  filter(Model == "Full model (243)") %>%
  normalize_categories()
cv_predictions <- read_csv(
  "results/tables/ten_fold_cv_predictions.csv",
  show_col_types = FALSE
) %>%
  filter(Model == "Full model (243)") %>%
  normalize_categories()

stopifnot(
  nrow(training_source) == 243,
  nrow(training_predictions) == 243,
  nrow(cv_predictions) == 243
)

predictive_metrics <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  tibble(
    n = length(observed),
    Predictive_R2 = 1 - sum((observed - predicted)^2) /
      sum((observed - mean(observed))^2),
    Pearson_r2 = cor(observed, predicted)^2,
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted)),
    Bias = mean(predicted - observed)
  )
}

training_metric <- predictive_metrics(
  training_predictions$Observed_LRV,
  training_predictions$Predicted_LRV
) %>%
  mutate(Evaluation = "Apparent training fit", .before = 1)

cv_metric <- predictive_metrics(
  cv_predictions$Observed_LRV,
  cv_predictions$Predicted_LRV
) %>%
  mutate(
    Evaluation = "Pooled out-of-fold predictions",
    .before = 1
  )

overall_metrics <- bind_rows(training_metric, cv_metric)

fold_metrics <- cv_predictions %>%
  group_by(Fold) %>%
  group_modify(
    ~ predictive_metrics(
      .x$Observed_LRV,
      .x$Predicted_LRV
    )
  ) %>%
  ungroup()

fold_composition <- cv_predictions %>%
  count(Fold, Disinfection_Type, name = "n") %>%
  complete(
    Fold = seq_len(10),
    Disinfection_Type = factor(
      DISINFECTION_LEVELS,
      levels = DISINFECTION_LEVELS
    ),
    fill = list(n = 0)
  ) %>%
  arrange(Fold, Disinfection_Type)

source_data <- cv_predictions %>%
  select(
    Fold,
    Site,
    Pathogen,
    Disinfection_Type,
    Observed_LRV,
    Predicted_LRV
  ) %>%
  left_join(
    training_source %>%
      transmute(
        Site,
        Pathogen,
        Disinfection_Type,
        Delta_Breadth,
        Delta_SNV_Density,
        Delta_Pi,
        Initial_Breadth = breadth_inf
      ),
    by = c("Site", "Pathogen", "Disinfection_Type")
  ) %>%
  mutate(
    Prediction_Error = Predicted_LRV - Observed_LRV,
    Absolute_Error = abs(Prediction_Error),
    Squared_Error = Prediction_Error^2
  ) %>%
  arrange(Fold, Disinfection_Type, Site, Pathogen)

stopifnot(
  nrow(source_data) == 243,
  all(
    complete.cases(
      source_data %>%
        select(
          Delta_Breadth,
          Delta_SNV_Density,
          Delta_Pi,
          Initial_Breadth
        )
    )
  )
)

write_csv(
  source_data,
  file.path(
    output_dir,
    "data",
    "Full_Model_243_TenFold_Prediction_Source_Data.csv"
  )
)
write_csv(
  overall_metrics,
  file.path(
    output_dir,
    "data",
    "Full_Model_243_Overall_Metrics.csv"
  )
)
write_csv(
  fold_metrics,
  file.path(
    output_dir,
    "data",
    "Full_Model_243_Fold_Metrics.csv"
  )
)
write_csv(
  fold_composition,
  file.path(
    output_dir,
    "data",
    "Full_Model_243_Fold_Composition.csv"
  )
)

cv_label <- sprintf(
  paste0(
    "n = %d\n",
    "Predictive R² = %.3f\n",
    "Pearson r² = %.3f\n",
    "RMSE = %.3f, MAE = %.3f\n",
    "Bias = %.3f"
  ),
  cv_metric$n,
  cv_metric$Predictive_R2,
  cv_metric$Pearson_r2,
  cv_metric$RMSE,
  cv_metric$MAE,
  cv_metric$Bias
)

p_cv <- ggplot(
  cv_predictions,
  aes(Observed_LRV, Predicted_LRV)
) +
  geom_point(
    aes(
      color = Disinfection_Type,
      size = breadth_inf
    ),
    alpha = 0.72
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  geom_smooth(
    method = "lm",
    color = "black",
    fill = "grey75",
    alpha = 0.30,
    linewidth = 0.8
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = cv_label,
    hjust = -0.05,
    vjust = 1.04,
    size = 3.7,
    lineheight = 0.96,
    family = MODEL_FONT,
    fontface = "bold",
    fill = "white",
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = ORIGINAL_COLORS,
    drop = FALSE
  ) +
  scale_size_continuous(range = c(2, 7)) +
  theme_original_model() +
  theme(
    aspect.ratio = 1,
    legend.position = "right"
  ) +
  labs(
    x = "Observation (True genomic LRV)",
    y = "Prediction (Predicted genomic LRV)",
    color = "Disinfection Process",
    size = "Initial Breadth"
  )

figure_base <- file.path(
  output_dir,
  "figure",
  "Fig_Sx_TenFold_CV_Full_Model_243"
)
save_original_figure(
  p_cv,
  figure_base,
  width = 7.5,
  height = 6,
  dpi = 600
)

print(overall_metrics)
print(fold_metrics)
