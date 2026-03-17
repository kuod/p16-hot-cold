# r_packages.R
# One-time R package installation for the cold-tumor-p16 analysis.
# Run this once before knitting the .Rmd notebooks:
#   Rscript r_packages.R

pkgs <- c(
  # Data I/O
  "arrow",      # read parquet files from Python notebooks
  "jsonlite",   # fromJSON() for tcga_tmb.json (notebook 05)
  # Data wrangling
  "dplyr",      # data manipulation
  "tidyr",      # pivot_wider for heatmap matrix
  "tibble",     # rownames_to_column / column_to_rownames
  "purrr",      # map_dfr for tidy iteration
  # Visualisation
  "ggplot2",    # plotting
  "pheatmap",   # Spearman correlation heatmap
  # Survival analysis
  "survival",   # Surv(), survfit(), survdiff(), coxph(), cox.zph()
  "survminer",  # ggsurvplot(), ggforest()
  "cmprsk",     # Fine-Gray competing risks regression (crr())
  # Statistical testing
  "broom",      # tidy() Cox model output
  "rstatix",    # cor_test() with BH FDR correction
  "car"         # vif() multicollinearity check for Cox models
  # Reproducibility
  "renv"        # lock-file based reproducibility (run renv::snapshot() after install)
)

failed <- character(0)
for (pkg in pkgs) {
  tryCatch(
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE),
    error = function(e) {
      failed <<- c(failed, pkg)
      message(sprintf("ERROR: Failed to install %s: %s", pkg, e$message))
    }
  )
}
if (length(failed) > 0) stop(sprintf("Install failed for: %s", paste(failed, collapse = ", ")))

# Verify all packages load cleanly and log versions
pkg_versions <- lapply(pkgs, function(p) {
  library(p, character.only = TRUE)
  v <- as.character(packageVersion(p))
  message(sprintf("OK: %-15s %s", p, v))
  data.frame(package = p, version = v, stringsAsFactors = FALSE)
})

version_df <- do.call(rbind, pkg_versions)
write.csv(version_df, "r_package_versions.csv", row.names = FALSE)
message("\nPackage versions saved to r_package_versions.csv")
message("Run renv::snapshot() to create renv.lock for reproducibility.")

# Auto-snapshot after install to lock package versions
if (requireNamespace("renv", quietly = TRUE)) {
  message("Running renv::snapshot() to lock package versions...")
  renv::snapshot(prompt = FALSE)
  message("renv.lock written.")
}
