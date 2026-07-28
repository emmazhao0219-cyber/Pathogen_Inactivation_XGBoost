rm(list = ls())

# Make the workflow independent of the caller's current working directory.
command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  # Rscript may encode spaces as "~+~" in commandArgs() on macOS.
  script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
  setwd(dirname(normalizePath(script_path)))
}

required_packages <- c(
  "readr", "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "xgboost", "ggplot2", "svglite", "ragg", "ggbeeswarm",
  "systemfonts"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Run `Rscript scripts/00_check_environment.R` for details."
  )
}

dir.create("models", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

source("scripts/01_model_training.R")
source("scripts/02_blind_validation_strategies.R")
source("scripts/03_data_visualization.R")
source("scripts/04_shap_full_model.R")
source("scripts/05_prediction_error_by_disinfection.R")

if (identical(Sys.getenv("BUILD_SI"), "1")) {
  source("scripts/06_si_internal_validation_material.R")
}

capture.output(
  sessionInfo(),
  file = "results/sessionInfo.txt"
)

cat("\nThree-strategy XGBoost validation completed successfully.\n")
