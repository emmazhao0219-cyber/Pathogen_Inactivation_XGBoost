# ==============================================================================
# Script 03: Academic Data Visualization
# ==============================================================================
rm(list = ls())

library(tidyverse)
library(readxl)
library(ggpubr)
library(xgboost)
library(extrafont)


if(!("Arial" %in% fonts())) {
  font_import(prompt = FALSE)
  loadfonts(device = "all")
}
my_colors <- c("UV" = "#ffb4ac", "Cl" = "#679186", "O3" = "#264e70", 
               "UV_Cl" = "#e4e1ca", "O3_Cl" = "#FF7F00", "UF" = "#6A3D9A")

cat("\n📊 Starting Visualization Pipeline...\n")

# ------------------------------------------------------------------------------
# 1. Feature Gain Barplot
# ------------------------------------------------------------------------------
cat("[1/4] Generating Feature Gain Plot...\n")
final_model <- xgb.load("models/XGBoost_Species_Agnostic_v1.model")
dv_encoder <- readRDS("models/Feature_Encoder_dv.rds")

dummy_data <- data.frame(LRV=0, Delta_Breadth=0, Delta_SNV_Density=0, Delta_Pi=0, Disinfection_Type="UV", breadth_inf=0)
dummy_matrix <- predict(dv_encoder, newdata = dummy_data)

imp_matrix <- xgb.importance(feature_names = colnames(dummy_matrix), model = final_model)
imp_df <- as.data.frame(imp_matrix) %>%
  mutate(
    Feature_Category = case_when(
      grepl("Disinfection_Type", Feature) ~ "Disinfection Context (Process Type)",
      Feature == "Relative_Breadth_Loss" ~ "Physical Damage (Relative Breadth Loss)",
      Feature == "Delta_Breadth" ~ "Physical Damage (Delta Breadth)",
      Feature == "Combined_Genetic_Damage" ~ "Genetic Damage (Combined SNV & Pi)",
      Feature == "Delta_SNV_Density" ~ "Genetic Damage (Delta SNV Density)",
      Feature == "Delta_Pi" ~ "Genetic Damage (Delta Pi)",
      Feature == "breadth_inf" ~ "Initial State (Breadth Inf)",
      TRUE ~ "Other Features"
    )
  ) %>%
  group_by(Feature_Category) %>% summarize(Total_Gain = sum(Gain), .groups = 'drop') %>%
  mutate(Contribution_Pct = Total_Gain * 100) %>% arrange(Contribution_Pct) %>%
  mutate(Feature_Category = factor(Feature_Category, levels = Feature_Category))

p_importance <- ggplot(imp_df, aes(x = Feature_Category, y = Contribution_Pct)) +
  geom_bar(stat = "identity", width = 0.6, fill = "#2C3E50") +
  geom_text(aes(label = sprintf("%.1f%%", Contribution_Pct)), hjust = -0.2, size = 4.5, fontface = "bold", family = "Arial", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(imp_df$Contribution_Pct) * 1.15), expand = c(0, 0)) +
  theme_bw() + theme(text = element_text(family = "Arial", face = "bold"), axis.title = element_text(size = 14), axis.text = element_text(size = 12), panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x = "Predictors", y = "Relative Feature Importance (Gain %)")

ggsave("results/figures/Feature_Importance_Gain.pdf", plot = p_importance, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 2. Training Set Fit Multi-Mapped
# ------------------------------------------------------------------------------
cat("[2/4] Generating Training Set Fit Scatter Plot...\n")
train_data <- read.csv("data/processed/training_data_ready.csv")
x_train <- predict(dv_encoder, newdata = train_data)
dtrain <- xgb.DMatrix(data = as.matrix(x_train), label = train_data$LRV)
train_data$Prediction <- predict(final_model, dtrain)

rmse_train <- sqrt(mean((train_data$Prediction - train_data$LRV)^2))

p_train_multi <- ggplot(train_data, aes(x = LRV, y = Prediction)) +
  geom_point(aes(color = Disinfection_Type, size = breadth_inf), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey", alpha = 0.8, linewidth = 0.8) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.x.npc = "left", label.y.npc = 0.90, size = 5, family = "Arial", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = sprintf("Train~RMSE == %.4f", rmse_train), parse = TRUE, hjust = -0.15, vjust = 8.5, size = 5, family = "Arial", fontface = "bold") +
  scale_color_manual(values = my_colors) + scale_size_continuous(range = c(2, 7)) +
  theme_bw() + theme(text = element_text(family = "Arial", face = "bold"), aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), legend.position = "right") +
  labs(x = "Observation (True LRV)", y = "Prediction (Predicted LRV)", color = "Disinfection Process", size = "Initial Breadth")

ggsave("results/figures/Training_Set_Fit_MultiMapped.pdf", plot = p_train_multi, width = 7.5, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Faceted Prediction Plot
# ------------------------------------------------------------------------------
cat("[3/4] Generating Faceted Scatter Plot...\n")
plot_df <- read_excel("results/tables/LRV_True_vs_Predicted.xlsx")

metrics_df <- plot_df %>% group_by(Disinfection_Type) %>% summarize(R2 = cor(LRV, Predicted_LRV)^2, RMSE = sqrt(mean((Predicted_LRV - LRV)^2)), .groups = 'drop') %>% mutate(rmse_label = sprintf("RMSE == %.3f", RMSE))

p_facet <- ggplot(plot_df, aes(x = LRV, y = Predicted_LRV)) +
  geom_point(aes(color = Disinfection_Type), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey80", alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~ Disinfection_Type, scales = "fixed") + 
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.x.npc = "left", label.y.npc = "top", family = "Arial", fontface = "bold", size = 4) +
  geom_text(data = metrics_df, aes(x = -Inf, y = Inf, label = rmse_label), hjust = -0.09, vjust = 4.3, parse = TRUE, family = "Arial", fontface = "bold", size = 4) +
  scale_color_manual(values = my_colors) + theme_bw() + theme(text = element_text(family = "Arial", face = "bold"), aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1.2), strip.text = element_text(size = 13, face = "bold"), legend.position = "none") +
  labs(x = "Observation (True LRV)", y = "Prediction (Predicted LRV)")

ggsave("results/figures/Prediction_Faceted_by_Disinfection.pdf", plot = p_facet, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. Absolute Error Boxplot
# ------------------------------------------------------------------------------
cat("[4/4] Generating Absolute Error Diagnosis Boxplot...\n")
error_df <- plot_df %>% mutate(Abs_Error = abs(Predicted_LRV - LRV))
summary_df <- error_df %>% group_by(Disinfection_Type) %>% summarize(MAE = mean(Abs_Error)) %>% arrange(MAE)
error_df$Disinfection_Type <- factor(error_df$Disinfection_Type, levels = summary_df$Disinfection_Type)

p_error <- ggplot(error_df, aes(x = Disinfection_Type, y = Abs_Error, fill = Disinfection_Type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black", linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5, color = "grey20", stroke = 0.5) +
  geom_text(data = summary_df, aes(x = Disinfection_Type, y = max(error_df$Abs_Error) * 1.08, label = sprintf("MAE\n%.2f", MAE)), size = 4.5, fontface = "bold", family = "Arial", color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = my_colors) + theme_bw() + theme(text = element_text(family = "Arial", face = "bold"), axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13), panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), legend.position = "none") +
  labs(x = "Disinfection Process (Arranged by Error Margin)", y = "Absolute Prediction Error (Log units)") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggsave("results/figures/Absolute_Error_Diagnosis.pdf", plot = p_error, width = 8, height = 6, dpi = 300)

cat("\n✅ All visualizations successfully generated in 'results/figures/'!\n")