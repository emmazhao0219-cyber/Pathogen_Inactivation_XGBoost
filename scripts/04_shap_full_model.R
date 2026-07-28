rm(list = ls())
source("scripts/model_utils.R")
suppressPackageStartupMessages(library(ggbeeswarm))

# Figure contract:
# The 243-record full model is primarily driven by Delta Breadth, while the
# direction and spread of every predictor's TreeSHAP contribution remain visible.
# This is a single-panel quantitative explanation figure. All contributions are
# calculated on the same complete-case training records used to fit the model.

training_full <- read_csv(
  "data/processed/training_full_model_completecase_243.csv",
  show_col_types = FALSE
) %>%
  normalize_categories()
assert_complete_features(
  training_full,
  FULL_FEATURES,
  "SHAP training data"
)

model_full <- xgb.load(
  "models/XGBoost_full_completecase_243.ubj"
)
feature_matrix <- make_feature_matrix(
  training_full,
  FULL_FEATURES
)
shap_raw <- predict(
  model_full,
  xgb.DMatrix(feature_matrix),
  predcontrib = TRUE
)

stopifnot(
  nrow(shap_raw) == nrow(training_full),
  ncol(shap_raw) == ncol(feature_matrix) + 1L
)
colnames(shap_raw) <- c(colnames(feature_matrix), "BIAS")

disinfection_columns <- grep(
  "^Disinfection_Type[.]",
  colnames(feature_matrix),
  value = TRUE
)

shap_grouped <- tibble(
  .row_id = seq_len(nrow(training_full)),
  `Delta Breadth` = shap_raw[, "Delta_Breadth"],
  `Disinfection Type` = rowSums(
    shap_raw[, disinfection_columns, drop = FALSE]
  ),
  `Delta Pi` = shap_raw[, "Delta_Pi"],
  `Delta SNV Density` = shap_raw[, "Delta_SNV_Density"],
  `Initial Breadth` = shap_raw[, "breadth_inf"]
)

value_grouped <- tibble(
  .row_id = seq_len(nrow(training_full)),
  `Delta Breadth` = training_full$Delta_Breadth,
  `Disinfection Type` = NA_real_,
  `Delta Pi` = training_full$Delta_Pi,
  `Delta SNV Density` = training_full$Delta_SNV_Density,
  `Initial Breadth` = training_full$breadth_inf
)

shap_long <- shap_grouped %>%
  pivot_longer(
    cols = -.row_id,
    names_to = "Feature",
    values_to = "SHAP_value"
  ) %>%
  left_join(
    value_grouped %>%
      pivot_longer(
        cols = -.row_id,
        names_to = "Feature",
        values_to = "Feature_value"
      ),
    by = c(".row_id", "Feature")
  ) %>%
  group_by(Feature) %>%
  mutate(
    Feature_value_scaled = if (
      all(is.na(Feature_value)) ||
      diff(range(Feature_value, na.rm = TRUE)) == 0
    ) {
      NA_real_
    } else {
      (Feature_value - min(Feature_value, na.rm = TRUE)) /
        diff(range(Feature_value, na.rm = TRUE))
    }
  ) %>%
  ungroup()

importance_summary <- shap_long %>%
  group_by(Feature) %>%
  summarise(
    mean_abs_SHAP = mean(abs(SHAP_value)),
    median_abs_SHAP = median(abs(SHAP_value)),
    min_SHAP = min(SHAP_value),
    max_SHAP = max(SHAP_value),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_SHAP))

feature_order <- importance_summary$Feature
shap_long <- shap_long %>%
  mutate(
    Feature = factor(
      Feature,
      levels = rev(feature_order)
    )
  )

write_csv(
  shap_long %>%
    mutate(Feature = as.character(Feature)) %>%
    arrange(.row_id, Feature),
  "results/tables/SHAP_full_model_243_source_data.csv"
)
write_csv(
  importance_summary,
  "results/tables/SHAP_full_model_243_importance.csv"
)

set.seed(2026)
p_shap <- ggplot(
  shap_long,
  aes(x = Feature, y = SHAP_value, color = Feature_value_scaled)
) +
  ggbeeswarm::geom_quasirandom(
    groupOnX = TRUE,
    width = 0.36,
    size = 1.55,
    alpha = 0.88,
    stroke = 0
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey45",
    linewidth = 0.45
  ) +
  coord_flip() +
  scale_color_gradientn(
    colors = c("#0066FF", "#A000B5", "#FF004F"),
    values = c(0, 0.5, 1),
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Low", "High"),
    na.value = "grey55",
    guide = guide_colorbar(
      title.position = "right",
      title.hjust = 0.5,
      barheight = grid::unit(42, "mm"),
      barwidth = grid::unit(2.5, "mm"),
      ticks = FALSE,
      frame.colour = NA
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  theme_classic(base_family = MODEL_FONT, base_size = 10) +
  theme(
    text = element_text(family = MODEL_FONT, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.65
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 9,
      angle = 90,
      vjust = 0.5
    ),
    legend.text = element_text(size = 8),
    plot.tag = element_text(
      family = MODEL_FONT,
      face = "bold",
      size = 12,
      hjust = 1.2,
      vjust = 0
    ),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(9, 7, 7, 13)
  ) +
  labs(
    tag = "e",
    x = NULL,
    y = "SHAP value (impact on Predicted genomic LRV)",
    color = "Feature value"
  )

save_original_figure(
  p_shap,
  "results/figures/SHAP_Summary_Full_Model_243",
  width = 7.2,
  height = 4.2,
  dpi = 600
)

print(importance_summary)
