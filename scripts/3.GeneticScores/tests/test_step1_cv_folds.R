#' ---------------------------
#'
#' Purpose: Tests for CV fold creation (Step 1)
#'
#' This script tests the CV fold creation functionality to ensure
#' it works correctly before proceeding to Step 2.
#'
#' ---------------------------
#'
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries
library(tidyverse)
library(survival)
library(caret)
library(testthat)

# Source the CV functions
source("scripts/3.GeneticScores/utils/create_cv_folds.R")

#' Test 1: Function loads and basic functionality works
test_basic_functionality <- function() {
  cat("=== TEST 1: Basic Functionality ===\n")
  
  # Create test data
  set.seed(42)
  test_data <- data.frame(
    ID = paste0("P", 1:100),
    LFS_YEARS = runif(100, 0.5, 10),
    LFS_STATUS = rbinom(100, 1, 0.3),
    AGE = runif(100, 20, 80),
    SEX = sample(c("M", "F"), 100, replace = TRUE)
  )
  
  # Test the function
  result <- create_cv_folds(test_data, k = 5, seed = 27)
  
  # Basic checks
  expect_equal(nrow(result), nrow(test_data))
  expect_true("fold" %in% colnames(result))
  expect_true(all(!is.na(result$fold)))
  expect_equal(length(unique(result$fold)), 5)
  expect_true(all(result$fold %in% 0:4))  # 0-based indexing
  
  cat("✓ Basic functionality works\n")
  cat("✓ All patients assigned to folds\n")
  cat("✓ Correct number of folds created\n")
  cat("✓ 0-based indexing working\n\n")
  
  return(result)
}

#' Test 2: Reproducibility with same seed
test_reproducibility <- function() {
  cat("=== TEST 2: Reproducibility ===\n")
  
  # Create test data
  test_data <- data.frame(
    ID = paste0("P", 1:50),
    LFS_YEARS = runif(50, 0.5, 10),
    LFS_STATUS = rbinom(50, 1, 0.3)
  )
  
  # Create folds twice with same seed
  result1 <- create_cv_folds(test_data, k = 5, seed = 27)
  result2 <- create_cv_folds(test_data, k = 5, seed = 27)
  
  # Check if results are identical
  expect_equal(result1$fold, result2$fold)
  
  cat("✓ Same seed produces identical results\n\n")
}

#' Test 3: Balance validation
test_balance <- function() {
  cat("=== TEST 3: Fold Balance ===\n")
  
  # Create test data with realistic medical distribution (similar to MDS data)
  set.seed(42)
  test_data <- data.frame(
    ID = paste0("P", 1:200),
    LFS_YEARS = rexp(200, rate = 0.3) + 0.1,  # Exponential survival times
    LFS_STATUS = rbinom(200, 1, 0.4)  # ~40% event rate (typical for MDS)
  )
  
  # Create folds
  result <- create_cv_folds(test_data, k = 10, seed = 27)
  
  # Validate balance
  validation <- validate_cv_folds(result)
  
  # Check balance metrics - más realista para datos médicos
  expect_true(validation$size_balance <= 10)  # Hasta 10 pacientes de diferencia
  expect_true(validation$event_rate_range <= 0.4)  # Hasta 40% de variación en eventos
  
  cat("✓ Fold sizes are reasonably balanced\n")
  cat("✓ Event rates are reasonably balanced\n")
  cat("Size balance:", validation$size_balance, "\n")
  cat("Event rate range:", round(validation$event_rate_range, 3), "\n\n")
  
  return(validation)
}

#' Test 4: Real data compatibility
test_real_data_compatibility <- function() {
  cat("=== TEST 4: Real Data Compatibility ===\n")
  
  # Check if patient_variables.tsv exists
  input_file <- "results/gnn/preprocess/patient_variables.tsv"
  
  if (!file.exists(input_file)) {
    cat("⚠ Real data file not found, skipping test\n")
    cat("Expected file:", input_file, "\n\n")
    return(NULL)
  }
  
  # Load real data
  real_data <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  cat("Real data loaded:", nrow(real_data), "patients\n")
  
  # Test with real data
  result <- create_cv_folds(real_data, k = 10, seed = 27)
  validation <- validate_cv_folds(result)
  
  # Basic checks
  expect_true(nrow(result) > 0)
  expect_equal(length(unique(result$fold)), 10)
  
  cat("✓ Real data processed successfully\n")
  cat("✓ 10 folds created\n")
  cat("Balance status:", validation$is_balanced, "\n\n")
  
  return(list(data = result, validation = validation))
}


#' Test 5: Edge cases
test_edge_cases <- function() {
  cat("=== TEST 5: Edge Cases ===\n")
  
  # Test with missing values
  test_data_na <- data.frame(
    ID = paste0("P", 1:20),
    LFS_YEARS = c(runif(15, 0.5, 10), rep(NA, 5)),
    LFS_STATUS = c(rbinom(15, 1, 0.3), rep(NA, 5))
  )
  
  result_na <- create_cv_folds(test_data_na, k = 3, seed = 27)
  expect_equal(nrow(result_na), 15)  # Should exclude NA rows
  
  # Test with very small data
  test_data_small <- data.frame(
    ID = paste0("P", 1:8),
    LFS_YEARS = runif(8, 0.5, 10),
    LFS_STATUS = rbinom(8, 1, 0.5)
  )
  
  result_small <- create_cv_folds(test_data_small, k = 3, seed = 27)
  expect_equal(length(unique(result_small$fold)), 3)
  
  cat("✓ Handles missing values correctly\n")
  cat("✓ Works with small datasets\n\n")
}

#' Run all tests
run_all_tests <- function() {
  cat("##########################################\n")
  cat("# TESTING CV FOLD CREATION (STEP 1)     #\n")
  cat("##########################################\n\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    test_basic_functionality()
    test_reproducibility()
    test_balance()
    real_data_result <- test_real_data_compatibility()
    test_file_operations()
    test_edge_cases()
    
    end_time <- Sys.time()
    
    cat("##########################################\n")
    cat("# ALL TESTS PASSED ✓                    #\n")
    cat("##########################################\n")
    cat("Test duration:", round(as.numeric(end_time - start_time), 2), "seconds\n\n")
    
    
  }, error = function(e) {
    cat("##########################################\n")
    cat("# TEST FAILED ✗                         #\n")
    cat("##########################################\n")
    cat("Error:", e$message, "\n")
    stop(e)
  })
}

#' Quick test for development
test_quick <- function() {
  cat("=== QUICK TEST ===\n")
  
  # Just test basic functionality
  test_data <- data.frame(
    ID = paste0("P", 1:30),
    LFS_YEARS = runif(30, 0.5, 10),
    LFS_STATUS = rbinom(30, 1, 0.3)
  )
  
  result <- create_cv_folds(test_data, k = 5, seed = 27)
  validation <- validate_cv_folds(result)
  
  cat("Quick test completed ✓\n")
  cat("Folds balanced:", validation$is_balanced, "\n")
  
  return(result)
}

# Uncomment to run all tests:
# run_all_tests()

# Uncomment for quick test:
# test_quick()
