#' ---------------------------
#'
#' Purpose: Execute and validate Step 1 - CV Fold Creation
#'
#' This script runs the CV fold creation and validates it works correctly
#' with your actual data before proceeding to Step 2.
#'
#' ---------------------------
#'
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

cat("##########################################\n")
cat("# STEP 1: CV FOLD CREATION               #\n")
cat("##########################################\n\n")

# Load libraries
library(tidyverse)
library(survival)
library(caret)

# Source functions
source("scripts/3.GeneticScores/utils/create_cv_folds.R")

load("results/preprocess/clinical_preproc.Rdata")

## Filter clinical data
cat("\n=== Filtering clinical data ===\n")
# Exclude patients with missing survival data and survival equal to 0
clinical <- clinical %>%
  filter(!is.na(LFS_YEARS) & !is.na(LFS_STATUS)) %>%
  filter(LFS_YEARS > 0)
# Create CV folds
cat("\n=== Creating 10-fold CV partitions ===\n")
clinical_fold <- create_cv_folds(clinical, time_col = "LFS_YEARS", 
  status_col = "LFS_STATUS", k = 10, seed = 27)

#' Step 1B: Validate fold quality
cat("\n=== Step 1B: Validating fold quality ===\n")
validation <- validate_cv_folds(clinical_fold, 
  time_col = "LFS_YEARS", status_col = "LFS_STATUS")

cat("\nValidation Results:\n")
cat("- Fold size balance (max - min):", validation$size_balance, "\n")
cat("- Event rate range (max - min):", round(validation$event_rate_range, 3), "\n")
cat("- Chi-square test p-value:", round(validation$chi_square_p, 3), "\n")
cat("- Folds are well balanced:", validation$is_balanced, "\n")

# Interpretación de resultados
if (validation$size_balance <= 5) {
  cat("✓ Fold sizes are well balanced (≤5 patients difference)\n")
} else if (validation$size_balance <= 10) {
  cat("~ Fold sizes are reasonably balanced (≤10 patients difference)\n")
} else {
  cat("⚠ Fold sizes show some imbalance (>10 patients difference)\n")
}

if (validation$event_rate_range <= 0.2) {
  cat("✓ Event rates are well balanced (≤20% variation)\n")
} else if (validation$event_rate_range <= 0.4) {
  cat("~ Event rates are reasonably balanced (≤40% variation)\n")
} else {
  cat("⚠ Event rates show significant variation (>40%)\n")
}

if (validation$is_balanced) {
  cat("✓ Overall fold quality is GOOD\n")
} else {
  cat("~ Overall fold quality is ACCEPTABLE for medical data\n")
}

#' Step 1C: Save CV folds for next steps
cat("\n=== Step 1C: Saving CV folds ===\n")

save(clinical_fold, file = "results/gnn/preprocess/clinical_fold.Rdata")
save(validation, file = "results/gnn/preprocess/clinical_fold_validation.Rdata")
write.table(clinical_fold, 
  file = "results/gnn/preprocess/clinical_fold.tsv", 
  sep = "\t", row.names = FALSE, quote = FALSE)

# Verify reproducibility
data_check <- create_cv_folds(clinical, k = 10, seed = 27)
if (identical(clinical_fold$fold, data_check$fold)) {
  cat("✓ CV fold creation is reproducible\n")
} else {
  stop("✗ CV fold creation is not reproducible")
}

#' Step 1D: Summary and next steps
cat("\n##########################################\n")
cat("# DEFINITION OF 10-FOLDS COMPLETED SUCCESSFULLY ✓       #\n")
cat("##########################################\n")

cat("\nSummary:\n")
cat("- Created 10-fold CV partitions from", nrow(clinical), "patients\n")
cat("- Excluded", nrow(clinical) - nrow(clinical_fold), "patients with missing survival data\n")
cat("- Final dataset:", nrow(clinical_fold), "patients with complete data\n")
cat("- Fold balance quality:", ifelse(validation$is_balanced, "GOOD", "ACCEPTABLE"), "\n")
