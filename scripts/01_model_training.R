# ==============================================================================
# Script 01: Data Preprocessing & XGBoost Engine Training
# ==============================================================================
rm(list = ls())

##loading package
library(readxl)
library(tidyverse)
library(xgboost)
library(caret)

cat("\n[1/4] Loading data...\n")

df <- read.csv("data/processed/training_data_ready.csv", stringsAsFactors = FALSE)

##Species-Agnostic Learning
model_ready_df <- df %>%
  select(LRV, Delta_Breadth, Delta_SNV_Density, Delta_Pi, Disinfection_Type, breadth_inf) %>%
  na.omit()

cat("\n[2/4] Performing feature engineering (One-hot Encoding)...\n")
#One-Hot encoding and DMatrix transfer
dv <- dummyVars(LRV ~ ., data = model_ready_df)
x_matrix <- predict(dv, newdata = model_ready_df)
y_label <- model_ready_df$LRV
dtrain <- xgb.DMatrix(data = as.matrix(x_matrix), label = y_label)

cat("\n[3/4] Training XGBoost engine...\n")
##training
params <- list(
  objective = "reg:squarederror",
  eta = 0.1,
  max_depth = 4,
  min_child_weight = 3,
  subsample = 0.8
)

set.seed(2026)
final_model <- xgb.train(params = params, data = dtrain, nrounds = 100)

cv_results <- xgb.cv(params = params, data = dtrain, nrounds = 100, nfold = 10, metrics = list("rmse", "mae"), verbose = FALSE)
predictions <- predict(final_model, dtrain)
r_squared <- cor(y_label, predictions)^2
rmse <- sqrt(mean((y_label - predictions)^2))

cat("\n==========================================\n")
cat(sprintf("📈 Train R²:   %.4f\n", r_squared))
cat(sprintf("📉 Train RMSE: %.4f\n", rmse))
cat(sprintf("📋 10-fold CV RMSE: %.4f\n", min(cv_results$evaluation_log$test_rmse_mean)))
cat("==========================================\n")

cat("\n[4/4] Encapsulating and saving model...\n")
#save
xgb.save(final_model, "models/XGBoost_Species_Agnostic_v1.model")
saveRDS(dv, "models/Feature_Encoder_dv.rds")
cat("✅ Training complete. Model saved to 'models/' directory.\n")