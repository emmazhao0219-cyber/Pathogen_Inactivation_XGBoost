required_packages <- c(
  "readr", "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "xgboost", "ggplot2", "svglite", "ragg", "ggbeeswarm",
  "systemfonts"
)

available <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

cat("R version:", R.version.string, "\n")
cat("Required packages:\n")
for (i in seq_along(required_packages)) {
  cat(
    sprintf(
      "  %-14s %s\n",
      required_packages[[i]],
      if (available[[i]]) {
        as.character(packageVersion(required_packages[[i]]))
      } else {
        "MISSING"
      }
    )
  )
}

if (!all(available)) {
  missing <- required_packages[!available]
  stop(
    "Install missing packages before running the pipeline: ",
    paste(missing, collapse = ", ")
  )
}

fonts <- unique(systemfonts::system_fonts()$family)
cat(
  "Plot font:",
  if (any(tolower(fonts) == "arial")) "Arial" else "sans (Arial fallback)",
  "\n"
)
cat("Environment check passed.\n")

