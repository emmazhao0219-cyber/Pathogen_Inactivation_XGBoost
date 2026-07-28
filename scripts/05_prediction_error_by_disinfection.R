rm(list = ls())
source("scripts/model_utils.R")

# Figure contract:
# Under the strict no-imputation blind validation, show the full distribution
# of absolute genomic-LRV prediction errors across disinfection types, together
# with group MAE and multiplicity-adjusted pairwise comparisons.

strict_blind <- read_csv(
  "results/tables/blind_primary_no_imputation_predictions.csv",
  show_col_types = FALSE
) %>%
  normalize_categories() %>%
  mutate(
    Absolute_Prediction_Error = abs(
      Predicted_LRV - Observed_LRV
    )
  )

stopifnot(
  nrow(strict_blind) == 157,
  all(is.finite(strict_blind$Absolute_Prediction_Error))
)

group_summary <- strict_blind %>%
  group_by(Disinfection_Type) %>%
  summarise(
    n = n(),
    MAE = mean(Absolute_Prediction_Error),
    Median_AE = median(Absolute_Prediction_Error),
    .groups = "drop"
  )

pair_grid <- combn(
  DISINFECTION_LEVELS,
  2,
  simplify = FALSE
)
pairwise_raw <- bind_rows(lapply(pair_grid, function(pair) {
  x <- strict_blind %>%
    filter(as.character(Disinfection_Type) == pair[[1]]) %>%
    pull(Absolute_Prediction_Error)
  y <- strict_blind %>%
    filter(as.character(Disinfection_Type) == pair[[2]]) %>%
    pull(Absolute_Prediction_Error)
  test <- suppressWarnings(
    wilcox.test(x, y, exact = FALSE, alternative = "two.sided")
  )
  tibble(
    group1 = pair[[1]],
    group2 = pair[[2]],
    n1 = length(x),
    n2 = length(y),
    p_raw = test$p.value
  )
}))

pairwise_results <- pairwise_raw %>%
  mutate(
    p_BH = p.adjust(p_raw, method = "BH"),
    significance = case_when(
      p_BH < 0.0001 ~ "****",
      p_BH < 0.001 ~ "***",
      p_BH < 0.01 ~ "**",
      p_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x1 = match(group1, DISINFECTION_LEVELS),
    x2 = match(group2, DISINFECTION_LEVELS),
    span = x2 - x1
  ) %>%
  arrange(span, x1) %>%
  mutate(
    y = 1.00 + (row_number() - 1) * 0.080
  )

write_csv(
  strict_blind %>%
    mutate(Disinfection_Type = as.character(Disinfection_Type)),
  "results/tables/Prediction_Error_Strict_No_Imputation_Source_Data.csv"
)
write_csv(
  group_summary %>%
    mutate(Disinfection_Type = as.character(Disinfection_Type)),
  "results/tables/Prediction_Error_Strict_No_Imputation_Group_Summary.csv"
)
write_csv(
  pairwise_results %>%
    select(group1, group2, n1, n2, p_raw, p_BH, significance),
  "results/tables/Prediction_Error_Strict_No_Imputation_Pairwise_Wilcoxon.csv"
)

bracket_height <- 0.028
p_error <- ggplot(
  strict_blind,
  aes(Disinfection_Type, Absolute_Prediction_Error)
) +
  geom_boxplot(
    aes(fill = Disinfection_Type),
    width = 0.64,
    outlier.shape = NA,
    linewidth = 0.65,
    color = "black",
    alpha = 0.9
  ) +
  geom_point(
    position = position_jitter(
      width = 0.14,
      height = 0,
      seed = 2026
    ),
    shape = 21,
    size = 1.7,
    stroke = 0.30,
    color = "grey25",
    fill = "grey45",
    alpha = 0.72
  ) +
  geom_segment(
    data = pairwise_results,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.38
  ) +
  geom_segment(
    data = pairwise_results,
    aes(
      x = x1,
      xend = x1,
      y = y,
      yend = y - bracket_height
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.38
  ) +
  geom_segment(
    data = pairwise_results,
    aes(
      x = x2,
      xend = x2,
      y = y,
      yend = y - bracket_height
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.38
  ) +
  geom_text(
    data = pairwise_results,
    aes(
      x = (x1 + x2) / 2,
      y = y + 0.018,
      label = significance
    ),
    inherit.aes = FALSE,
    family = MODEL_FONT,
    fontface = "bold",
    size = 2.6,
    vjust = 0
  ) +
  geom_text(
    data = group_summary,
    aes(
      x = Disinfection_Type,
      y = 0.90,
      label = sprintf("MAE\n%.2f", MAE)
    ),
    inherit.aes = FALSE,
    family = MODEL_FONT,
    fontface = "bold",
    size = 3.1,
    lineheight = 0.92
  ) +
  scale_fill_manual(
    values = ORIGINAL_COLORS,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      UV = "UV",
      Cl = "Cl",
      O3 = "O3",
      UV_Cl = "UV-Cl",
      O3_Cl = "O3-Cl",
      UF = "UF"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 2.22),
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0, 0.015))
  ) +
  theme_classic(base_family = MODEL_FONT, base_size = 10) +
  theme(
    text = element_text(family = MODEL_FONT, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.55),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.65
    ),
    legend.position = "none",
    plot.title = element_text(
      size = 11,
      hjust = 0.5,
      face = "bold"
    ),
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
    tag = "d",
    title = "Prediction Error Distribution by Disinfection Type",
    x = "Disinfection Type",
    y = "Absolute Prediction Error (genomic LRV)"
  )

save_original_figure(
  p_error,
  "results/figures/Prediction_Error_By_Disinfection_Strict_No_Imputation",
  width = 6.3,
  height = 5.3,
  dpi = 600
)

print(group_summary)
print(
  pairwise_results %>%
    select(group1, group2, p_raw, p_BH, significance)
)
