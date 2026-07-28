suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(xgboost)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(systemfonts)
})

set.seed(2026)

AVAILABLE_FONTS <- unique(systemfonts::system_fonts()$family)
MODEL_FONT <- if (any(tolower(AVAILABLE_FONTS) == "arial")) {
  "Arial"
} else {
  "sans"
}

DISINFECTION_LEVELS <- c("UV", "Cl", "O3", "UV_Cl", "O3_Cl", "UF")
FULL_FEATURES <- c(
  "Delta_Breadth", "Delta_SNV_Density", "Delta_Pi", "breadth_inf"
)
REDUCED_FEATURES <- c(
  "Delta_Breadth", "Delta_SNV_Density", "breadth_inf"
)
MODEL_PARAMS <- list(
  objective = "reg:squarederror",
  eta = 0.1,
  max_depth = 4,
  min_child_weight = 3,
  subsample = 0.8,
  nthread = 1
)
NROUNDS <- 100

# Exact palette and visual grammar retained from the original project.
ORIGINAL_COLORS <- c(
  UV = "#ffb4ac",
  Cl = "#679186",
  O3 = "#264e70",
  UV_Cl = "#e4e1ca",
  O3_Cl = "#FF7F00",
  UF = "#6A3D9A"
)
STRATEGY_COLORS <- c(
  "Strict no-imputation" = "#264e70",
  "Imputation-assisted" = "#679186",
  "Reduced model (no Delta_Pi)" = "#FF7F00"
)

normalize_categories <- function(dat) {
  dat %>%
    mutate(
      Disinfection_Type = factor(
        Disinfection_Type,
        levels = DISINFECTION_LEVELS
      ),
      Pathogen = as.character(Pathogen)
    )
}

required_columns <- function(features) {
  c(features, "Disinfection_Type")
}

assert_complete_features <- function(dat, features, label = "dataset") {
  required <- required_columns(features)
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop(label, " is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  incomplete <- !complete.cases(dat[, required, drop = FALSE])
  if (any(incomplete)) {
    stop(label, " contains ", sum(incomplete), " incomplete feature rows.")
  }
  invisible(TRUE)
}

make_feature_matrix <- function(dat, features) {
  dat <- normalize_categories(dat)
  assert_complete_features(dat, features, "feature matrix input")
  if (any(is.na(dat$Disinfection_Type))) {
    stop("Unknown disinfection category in feature matrix input.")
  }

  numeric_matrix <- as.matrix(dat[, features, drop = FALSE])
  storage.mode(numeric_matrix) <- "double"
  dummy <- vapply(
    DISINFECTION_LEVELS[-1],
    function(level) as.numeric(dat$Disinfection_Type == level),
    numeric(nrow(dat))
  )
  if (nrow(dat) == 1) {
    dummy <- matrix(dummy, nrow = 1)
  }
  colnames(dummy) <- paste0(
    "Disinfection_Type.",
    DISINFECTION_LEVELS[-1]
  )
  cbind(numeric_matrix, dummy)
}

fit_xgboost <- function(dat, features) {
  assert_complete_features(dat, features, "training data")
  xgb.train(
    params = MODEL_PARAMS,
    data = xgb.DMatrix(
      make_feature_matrix(dat, features),
      label = dat$LRV_model
    ),
    nrounds = NROUNDS,
    verbose = 0
  )
}

predict_xgboost <- function(model, dat, features) {
  assert_complete_features(dat, features, "prediction data")
  predict(model, xgb.DMatrix(make_feature_matrix(dat, features)))
}

metric_row <- function(
    observed,
    predicted,
    model,
    strategy,
    group = "Overall",
    sample_basis = "Available data") {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  tibble(
    Model = model,
    Strategy = strategy,
    Sample_basis = sample_basis,
    Group = group,
    n = length(observed),
    R2 = if (
      length(observed) >= 3 &&
      sd(observed) > 0 &&
      sd(predicted) > 0
    ) {
      cor(observed, predicted)^2
    } else {
      NA_real_
    },
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted)),
    Bias = mean(predicted - observed)
  )
}

metrics_by_disinfection <- function(
    dat,
    model,
    strategy,
    sample_basis = "Available data") {
  bind_rows(
    metric_row(
      dat$Observed_LRV,
      dat$Predicted_LRV,
      model,
      strategy,
      "Overall",
      sample_basis
    ),
    bind_rows(lapply(DISINFECTION_LEVELS, function(group) {
      d <- dat %>%
        filter(as.character(Disinfection_Type) == group)
      if (nrow(d) == 0) {
        return(NULL)
      }
      metric_row(
        d$Observed_LRV,
        d$Predicted_LRV,
        model,
        strategy,
        group,
        sample_basis
      )
    }))
  )
}

make_stratified_folds <- function(dat, k = 10, seed = 2026) {
  set.seed(seed)
  fold_id <- integer(nrow(dat))
  for (group in levels(droplevels(dat$Disinfection_Type))) {
    idx <- which(dat$Disinfection_Type == group)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold_id
}

cross_validate_xgboost <- function(
    dat,
    features,
    model_name,
    k = 10,
    seed = 2026) {
  fold_id <- make_stratified_folds(dat, k = k, seed = seed)
  bind_rows(lapply(seq_len(k), function(fold) {
    train_fold <- dat[fold_id != fold, , drop = FALSE]
    test_fold <- dat[fold_id == fold, , drop = FALSE]
    fold_model <- fit_xgboost(train_fold, features)
    test_fold %>%
      transmute(
        Model = model_name,
        Fold = fold,
        Site,
        Pathogen,
        Disinfection_Type,
        breadth_inf,
        Observed_LRV = LRV_model,
        Predicted_LRV = predict_xgboost(
          fold_model,
          test_fold,
          features
        )
      )
  }))
}

# Only training-set statistics are used. Missingness flags are audit variables
# and are deliberately excluded from both model matrices.
impute_blind_from_training <- function(blind, training, features) {
  blind <- normalize_categories(blind) %>%
    mutate(.blind_row_id = row_number())
  training <- normalize_categories(training)
  audit_parts <- vector("list", length(features))

  for (i in seq_along(features)) {
    feature <- features[[i]]
    missing_flag <- is.na(blind[[feature]])

    group_lookup <- training %>%
      transmute(
        Disinfection_Type = as.character(Disinfection_Type),
        Pathogen,
        value = .data[[feature]]
      ) %>%
      group_by(Disinfection_Type, Pathogen) %>%
      summarise(
        group_median = median(value),
        group_n = n(),
        .groups = "drop"
      )

    pathogen_lookup <- training %>%
      transmute(Pathogen, value = .data[[feature]]) %>%
      group_by(Pathogen) %>%
      summarise(
        pathogen_median = median(value),
        pathogen_n = n(),
        .groups = "drop"
      )

    overall_median <- median(training[[feature]])
    joined <- blind %>%
      mutate(
        Disinfection_Type_character = as.character(Disinfection_Type)
      ) %>%
      left_join(
        group_lookup,
        by = c(
          "Disinfection_Type_character" = "Disinfection_Type",
          "Pathogen"
        )
      ) %>%
      left_join(pathogen_lookup, by = "Pathogen")

    source <- case_when(
      !missing_flag ~ "Observed",
      is.finite(joined$group_median) ~
        "Disinfection_x_Pathogen_training_median",
      is.finite(joined$pathogen_median) ~
        "Pathogen_training_median",
      TRUE ~ "Overall_training_median"
    )
    imputed_value <- ifelse(
      !missing_flag,
      blind[[feature]],
      ifelse(
        is.finite(joined$group_median),
        joined$group_median,
        ifelse(
          is.finite(joined$pathogen_median),
          joined$pathogen_median,
          overall_median
        )
      )
    )

    blind[[paste0(feature, "_missing_original")]] <-
      as.integer(missing_flag)
    blind[[feature]] <- imputed_value
    audit_parts[[i]] <- tibble(
      .blind_row_id = blind$.blind_row_id,
      Site = blind$Site,
      Pathogen = blind$Pathogen,
      Disinfection_Type = as.character(blind$Disinfection_Type),
      Feature = feature,
      Was_missing = missing_flag,
      Imputation_source = source,
      Final_value = imputed_value
    )
  }

  assert_complete_features(blind, features, "imputed blind data")
  list(
    data = blind %>% select(-.blind_row_id),
    audit = bind_rows(audit_parts)
  )
}

theme_original_model <- function() {
  theme_bw() +
    theme(
      text = element_text(family = MODEL_FONT, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.2
      ),
      panel.grid.minor = element_blank()
    )
}

save_original_figure <- function(
    plot,
    filename,
    width = 8,
    height = 6,
    dpi = 600) {
  svglite::svglite(
    paste0(filename, ".svg"),
    width = width,
    height = height
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(filename, ".pdf"),
    width = width,
    height = height,
    family = MODEL_FONT
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(filename, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  print(plot)
  dev.off()

  if (identical(Sys.getenv("EXPORT_TIFF"), "1")) {
    ragg::agg_tiff(
      paste0(filename, ".tiff"),
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    print(plot)
    dev.off()
  }
}
