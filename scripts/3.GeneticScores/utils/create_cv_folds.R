#' ---------------------------
#'
#' Purpose of script: Create CV folds for cross-validation
#'
#' This script creates 10-fold cross-validation partitions that will be
#' used consistently across CoxPH and GNN models to ensure fair comparison.
#'
#' ---------------------------
#'
#' Notes:
#' - Uses the same seed (27) as the original train/test split
#' - Ensures stratification by survival status for balanced folds
#' - Generates reproducible partitions
#'
#' ---------------------------

library(tidyverse)
library(survival)
library(caret)

#' Create stratified k-fold cross-validation partitions
#' 
#' @param data Data frame with survival data
#' @param time_col Column name for survival time
#' @param status_col Column name for event status
#' @param id_col Column name for patient ID
#' @param k Number of folds (default: 10)
#' @param seed Random seed for reproducibility (default: 27)
#' @return Data frame with original data plus 'fold' column (0-based indexing)
create_cv_folds <- function(data, time_col = "LFS_YEARS", status_col = "LFS_STATUS", 
                           id_col = "ID", k = 10, seed = 27) {
  
  # Validate inputs
  required_cols <- c(time_col, status_col, id_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Remove rows with missing survival data
  complete_data <- data %>%
    filter(!is.na(.data[[time_col]]) & !is.na(.data[[status_col]]))
  
  cat("Original data:", nrow(data), "patients\n")
  cat("Complete survival data:", nrow(complete_data), "patients\n")

  # Stratify by time and status
  stratify_var <- paste(cut(complete_data[[time_col]], breaks = 8, labels = FALSE), 
    complete_data[[status_col]], sep = "_")

  # Set seed for reproducibility
  set.seed(seed)
  
  # Create stratified folds using survival data
  # This ensures balanced distribution of events across folds
  fold_indices <- createFolds(stratify_var, k = k, list = TRUE, returnTrain = FALSE)
  
  # Initialize fold column
  complete_data$fold <- NA
  
  # Assign fold numbers (0-based for consistency with Python)
  for (i in seq_along(fold_indices)) {
    fold_num <- i - 1  # 0-based indexing
    complete_data$fold[fold_indices[[i]]] <- fold_num
  }
  
  # Verify all patients have been assigned a fold
  if (any(is.na(complete_data$fold))) {
    stop("Error: Some patients were not assigned to any fold")
  }
  

  return(complete_data)
}

#' Validate CV folds quality
#' 
#' @param data Data frame with fold assignments
#' @param time_col Column name for survival time
#' @param status_col Column name for event status
#' @return List with validation results
validate_cv_folds <- function(data, time_col = "LFS_YEARS", status_col = "LFS_STATUS") {
  

  # Overall statistics
  cat("\nOverall statistics:\n")
  cat("Total patients:", nrow(data), "\n")
  cat("Total events:", sum(data[[status_col]]), "\n")
  cat("Overall event rate:", round(mean(data[[status_col]]), 3), "\n")
  cat("Median survival time:", round(median(data[[time_col]]), 2), "\n")
  

  fold_stats <- data %>%
    group_by(fold) %>%
    summarise(
      n_patients = n(),
      n_events = sum(.data[[status_col]]),
      event_rate = mean(.data[[status_col]]),
      median_time = median(.data[[time_col]]),
      .groups = 'drop'
    )
  
  # Check balance
  size_balance <- max(fold_stats$n_patients) - min(fold_stats$n_patients)
  event_rate_range <- max(fold_stats$event_rate) - min(fold_stats$event_rate)
  
  # Chi-square test for event distribution
  event_counts <- table(data$fold, data[[status_col]])
  chi_test <- chisq.test(event_counts)
  
  validation_results <- list(
    fold_stats = fold_stats,
    size_balance = size_balance,
    event_rate_range = event_rate_range,
    chi_square_p = chi_test$p.value,
    is_balanced = size_balance <= 5 & event_rate_range <= 0.3 & chi_test$p.value > 0.001
  )
  
  return(validation_results)
}



#' Create CV folds optimized for medical data with potentially unbalanced events
#' 
#' @param data Data frame with survival data
#' @param time_col Column name for survival time
#' @param status_col Column name for event status
#' @param id_col Column name for patient ID
#' @param k Number of folds (default: 10)
#' @param seed Random seed for reproducibility (default: 27)
#' @param max_attempts Maximum attempts to create balanced folds (default: 5)
#' @return Data frame with original data plus 'fold' column (0-based indexing)
create_cv_folds_robust <- function(data, time_col = "LFS_YEARS", status_col = "LFS_STATUS", 
                                  id_col = "ID", k = 10, seed = 27, max_attempts = 5) {
  
  best_result <- NULL
  best_balance <- Inf
  
  for (attempt in 1:max_attempts) {
    
    # Use different seed for each attempt but keep reproducible
    attempt_seed <- seed + attempt - 1
    
    # Create folds
    result <- create_cv_folds(data, time_col, status_col, id_col, k, attempt_seed)
    
    # Evaluate balance
    validation <- validate_cv_folds(result, time_col, status_col)
    
    # Calculate combined balance score (lower is better)
    balance_score <- validation$size_balance + (validation$event_rate_range * 10)
    
    if (balance_score < best_balance) {
      best_balance <- balance_score
      best_result <- result
    }
    
    # If we achieve good balance, stop early
    if (validation$is_balanced) {
      cat("Achieved good balance on attempt", attempt, "\n")
      break
    }
  }
  
  cat("Best balance achieved: size_diff =", 
      validate_cv_folds(best_result, time_col, status_col)$size_balance,
      ", event_rate_range =", 
      round(validate_cv_folds(best_result, time_col, status_col)$event_rate_range, 3), "\n")
  
  return(best_result)
}

# Example usage (commented out for safety):
# data_with_folds <- create_and_save_cv_folds(
#   input_path = "results/gnn/preprocess/patient_variables.tsv",
#   output_path = "results/gnn/preprocess/patient_variables_cv.tsv",
#   k = 10,
#   seed = 27
# )
