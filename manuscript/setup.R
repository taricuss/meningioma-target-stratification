# manuscript/setup.R — one-time bootstrap (run in your interactive R Console in VS Code)
# Fixes: "there is no package called 'bookdown'" and other missing packages.
# Run ONCE in the R Console (click the file, press Ctrl+Shift+S, or run:
#   source("manuscript/setup.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_pkgs <- c(
  "tidyverse", "here", "flextable", "officer",
  "scales", "stringr", "knitr", "bookdown", "purrr", "rmarkdown"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% rownames(installed.packages()))]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, Ncpus = max(1L, parallel::detectCores() - 1L))
} else {
  message("All ", length(required_pkgs), " required packages already installed.")
}

# Verify all load cleanly
failed <- character(0)
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) failed <- c(failed, p)
}
if (length(failed) == 0) {
  message("\nAll packages verified OK. You can now knit paper_manuscript.Rmd.")
} else {
  stop("FAILED to load: ", paste(failed, collapse = ", "))
}
