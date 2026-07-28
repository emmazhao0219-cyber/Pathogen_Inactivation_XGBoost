rm(list = ls())
source("scripts/model_utils.R")

training <- read_csv(
  "results/tables/training_fit_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
cv <- read_csv(
  "results/tables/ten_fold_cv_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
training_metrics <- read_csv(
  "results/tables/training_and_cv_metrics.csv",
  show_col_types = FALSE
)
blind_primary <- read_csv(
  "results/tables/blind_primary_no_imputation_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
blind_imputed <- read_csv(
  "results/tables/blind_sensitivity_imputed_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
blind_reduced <- read_csv(
  "results/tables/blind_reduced_model_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
metrics_available <- read_csv(
  "results/tables/blind_metrics_available_samples.csv",
  show_col_types = FALSE
)
metrics_common <- read_csv(
  "results/tables/blind_metrics_common_complete_subset.csv",
  show_col_types = FALSE
)
metrics_input_status <- read_csv(
  "results/tables/imputation_assisted_metrics_by_input_status.csv",
  show_col_types = FALSE
)
importance <- read_csv(
  "results/tables/feature_importance_raw.csv",
  show_col_types = FALSE
)

make_fit_plot <- function(dat, metric) {
  correlation_test <- cor.test(dat$Observed_LRV, dat$Predicted_LRV)
  p_text <- if (correlation_test$p.value < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", correlation_test$p.value)
  }
  metric_label <- sprintf(
    "n = %d\nPearson r² = %.2f, %s\nRMSE = %.4f, MAE = %.4f",
    metric$n,
    metric$R2,
    p_text,
    metric$RMSE,
    metric$MAE
  )

  ggplot(dat, aes(Observed_LRV, Predicted_LRV)) +
    geom_point(
      aes(color = Disinfection_Type, size = breadth_inf),
      alpha = 0.7
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
      fill = "grey",
      alpha = 0.35,
      linewidth = 0.8
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = metric_label,
      hjust = -0.05,
      vjust = 1.06,
      size = 4.0,
      lineheight = 1.0,
      family = MODEL_FONT,
      fontface = "bold",
      fill = "white",
      linewidth = 0.3
    ) +
    scale_color_manual(values = ORIGINAL_COLORS, drop = FALSE) +
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
}

fit_specs <- tribble(
  ~Model, ~Split, ~File,
  "Full model (243)", "Apparent training fit",
  "Training_Set_Fit_Full_Model_243",
  "Full model (243)", "Stratified 10-fold CV",
  "TenFold_CV_Fit_Full_Model_243",
  "Reduced model (430)", "Apparent training fit",
  "Training_Set_Fit_Reduced_Model_430",
  "Reduced model (430)", "Stratified 10-fold CV",
  "TenFold_CV_Fit_Reduced_Model_430"
)

for (i in seq_len(nrow(fit_specs))) {
  spec <- fit_specs[i, ]
  dat <- if (spec$Split == "Apparent training fit") {
    training %>% filter(Model == spec$Model)
  } else {
    cv %>% filter(Model == spec$Model)
  }
  metric <- training_metrics %>%
    filter(Model == spec$Model, Strategy == spec$Split)
  p <- make_fit_plot(dat, metric)
  save_original_figure(
    p,
    file.path("results/figures", spec$File),
    width = 7.5,
    height = 6
  )
}

make_blind_facet <- function(dat, metrics) {
  labels <- metrics %>%
    filter(Group != "Overall") %>%
    transmute(
      Disinfection_Type = factor(Group, levels = DISINFECTION_LEVELS),
      metric_label = sprintf(
        "n = %d\nPearson r² = %.3f\nRMSE = %.3f\nMAE = %.3f",
        n,
        R2,
        RMSE,
        MAE
      )
    )

  ggplot(dat, aes(Observed_LRV, Predicted_LRV)) +
    geom_point(
      aes(color = Disinfection_Type),
      size = 3,
      alpha = 0.8
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
      fill = "grey80",
      alpha = 0.3,
      linewidth = 0.8
    ) +
    facet_wrap(~ Disinfection_Type, scales = "fixed") +
    geom_text(
      data = labels,
      aes(x = -Inf, y = Inf, label = metric_label),
      hjust = -0.05,
      vjust = 1.04,
      family = MODEL_FONT,
      fontface = "bold",
      size = 2.9,
      lineheight = 0.94,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = ORIGINAL_COLORS, drop = FALSE) +
    theme_original_model() +
    theme(
      aspect.ratio = 1,
      strip.background = element_rect(
        fill = "grey90",
        color = "black",
        linewidth = 1.2
      ),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      x = "Observation (True genomic LRV)",
      y = "Prediction (Predicted genomic LRV)"
    )
}

blind_specs <- list(
  list(
    data = blind_primary,
    strategy = "Strict no-imputation",
    file = "Blind_Primary_No_Imputation_Faceted"
  ),
  list(
    data = blind_imputed,
    strategy = "Imputation-assisted",
    file = "Blind_Sensitivity_Imputation_Assisted_Faceted"
  ),
  list(
    data = blind_reduced,
    strategy = "Reduced model (no Delta_Pi)",
    file = "Blind_Reduced_Model_No_DeltaPi_Faceted"
  )
)

for (spec in blind_specs) {
  metric <- metrics_available %>%
    filter(Strategy == spec$strategy)
  p <- make_blind_facet(spec$data, metric)
  save_original_figure(
    p,
    file.path("results/figures", spec$file),
    width = 8,
    height = 6
  )
}

plot_strategy_metrics <- function(metric_data) {
  strategy_labels <- c(
    "Strict no-imputation" = "Strict\ncomplete",
    "Imputation-assisted" = "Imputed\nfull",
    "Reduced model (no Delta_Pi)" = "Reduced\nno Delta Pi"
  )
  plot_df <- metric_data %>%
    filter(Group == "Overall") %>%
    mutate(
      Strategy = factor(
        Strategy,
        levels = names(STRATEGY_COLORS)
      )
    ) %>%
    pivot_longer(
      cols = c(R2, RMSE, MAE),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("R2", "RMSE", "MAE"),
        labels = c("Pearson r²", "RMSE", "MAE")
      ),
      value_label = sprintf("%.3f\n(n=%d)", Value, n)
    )

  ggplot(plot_df, aes(Strategy, Value, fill = Strategy)) +
    geom_col(width = 0.64, color = "black", linewidth = 0.6) +
    geom_text(
      aes(label = value_label),
      vjust = -0.25,
      size = 3.7,
      family = MODEL_FONT,
      fontface = "bold",
      lineheight = 0.95
    ) +
    facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = strategy_labels) +
    scale_fill_manual(values = STRATEGY_COLORS, drop = FALSE) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.22))
    ) +
    theme_original_model() +
    theme(
      axis.text.x = element_text(size = 9),
      strip.background = element_rect(
        fill = "grey90",
        color = "black",
        linewidth = 1.2
      ),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      x = "Validation strategy",
      y = "Performance metric"
    )
}

p_available <- plot_strategy_metrics(metrics_available)
save_original_figure(
  p_available,
  "results/figures/Validation_Strategies_Available_Samples_Comparison",
  width = 9,
  height = 4.8
)

p_common <- plot_strategy_metrics(metrics_common)
save_original_figure(
  p_common,
  "results/figures/Validation_Strategies_Common_Subset_Comparison",
  width = 9,
  height = 4.8
)

status_plot_data <- metrics_input_status %>%
  filter(Group == "Overall") %>%
  mutate(
    Input_status = if_else(
      str_detect(Strategy, "Complete$"),
      "Complete input",
      "Imputed input"
    )
  ) %>%
  pivot_longer(
    cols = c(R2, RMSE, MAE),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("R2", "RMSE", "MAE"),
      labels = c("Pearson r²", "RMSE", "MAE")
    ),
    value_label = sprintf("%.3f\n(n=%d)", Value, n)
  )

p_status <- ggplot(
  status_plot_data,
  aes(Input_status, Value, fill = Input_status)
) +
  geom_col(width = 0.60, color = "black", linewidth = 0.6) +
  geom_text(
    aes(label = value_label),
    vjust = -0.25,
    size = 3.8,
    family = MODEL_FONT,
    fontface = "bold",
    lineheight = 0.95
  ) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c(
      "Complete input" = "#264e70",
      "Imputed input" = "#679186"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Complete input" = "Complete",
      "Imputed input" = "Imputed"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22))
  ) +
  theme_original_model() +
  theme(
    axis.text.x = element_text(size = 11),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 1.2
    ),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Input status within imputation-assisted validation",
    y = "Performance metric"
  )

save_original_figure(
  p_status,
  "results/figures/Imputation_Assisted_Complete_vs_Imputed",
  width = 8,
  height = 4.8
)

importance_grouped <- importance %>%
  mutate(
    Feature_Category = case_when(
      grepl("Disinfection_Type", Feature) ~
        "Disinfection Context (Process Type)",
      Feature == "Delta_Breadth" ~
        "Physical Damage (Delta Breadth)",
      Feature == "Delta_SNV_Density" ~
        "Genetic Damage (Delta SNV Density)",
      Feature == "Delta_Pi" ~
        "Genetic Damage (Delta Pi)",
      Feature == "breadth_inf" ~
        "Initial State (Breadth Inf)",
      TRUE ~ "Other Features"
    )
  ) %>%
  group_by(Model, Feature_Category) %>%
  summarise(Total_Gain = sum(Gain), .groups = "drop") %>%
  group_by(Model) %>%
  mutate(Contribution_Pct = Total_Gain / sum(Total_Gain) * 100) %>%
  ungroup()
write_csv(
  importance_grouped,
  "results/tables/feature_importance_grouped.csv"
)

p_importance <- ggplot(
  importance_grouped,
  aes(
    reorder(Feature_Category, Contribution_Pct),
    Contribution_Pct,
    fill = Model
  )
) +
  geom_col(width = 0.62, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.1f%%", Contribution_Pct)),
    hjust = -0.15,
    size = 3.8,
    family = MODEL_FONT,
    fontface = "bold"
  ) +
  coord_flip() +
  facet_wrap(~ Model, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Full model (243)" = "#2C3E50",
      "Reduced model (430)" = "#FF7F00"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20))
  ) +
  theme_original_model() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 1.2
    ),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Predictors",
    y = "Relative Feature Importance (Gain %)"
  )

save_original_figure(
  p_importance,
  "results/figures/Feature_Importance_Model_Comparison",
  width = 8,
  height = 7
)
