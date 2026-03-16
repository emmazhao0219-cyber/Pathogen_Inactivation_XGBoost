# ==============================================================================
# Script 02: Inference Workflow for Unknown Datasets
# ==============================================================================
rm(list = ls())

library(tidyverse)
library(xgboost)
library(caret)
library(readxl)
library(writexl)

cat("\n[1/4] Loading inference data and pre-trained engine...\n")
##loading data
x_matrix <- read_excel("data/raw/data_no_internal_standard.xlsx")
final_model <- xgb.load("models/XGBoost_Species_Agnostic_v1.model")
dv_encoder <- readRDS("models/Feature_Encoder_dv.rds")

cat("\n[2/4] Aligning features...\n")

prediction_data <- x_matrix %>%
  mutate(
    Relative_Breadth_Loss = Delta_Breadth / (breadth_inf + 1e-5), 
    Combined_Genetic_Damage = Delta_SNV_Density * Delta_Pi,
    LRV = 0  # Placeholder
  )

prediction_data$Disinfection_Type <- factor(
  prediction_data$Disinfection_Type, 
  levels = c("UV", "Cl", "O3", "UV_Cl", "O3_Cl", "UF") 
)

cat("\n[3/4] Running XGBoost prediction...\n")
x_new_matrix <- predict(dv_encoder, newdata = prediction_data)
dnew <- xgb.DMatrix(data = as.matrix(x_new_matrix))
x_matrix$Predicted_LRV <- predict(final_model, dnew)

#print predicted results
final_output <- x_matrix %>%
  select(Site, Pathogen, Disinfection_Type, breadth_inf, Delta_Breadth, Predicted_LRV)
write_xlsx(final_output, "results/tables/LRV_Predicted_Output.xlsx")

cat("\n[4/4] Matching with True Observations...\n")
##comparision with true LRV data
df2_true <- read_excel("data/raw/original_true_LRV.xlsx") %>% select(Site, Pathogen, Disinfection_Type, LRV)
merged_df <- inner_join(df2_true, final_output %>% select(Site, Pathogen, Predicted_LRV), by = c("Site", "Pathogen")) %>%
  select(Site, Pathogen, Disinfection_Type, LRV, Predicted_LRV)

write_xlsx(merged_df, "results/tables/LRV_True_vs_Predicted.xlsx")
cat("✅ Inference complete. Result tables saved to 'results/tables/'.\n")