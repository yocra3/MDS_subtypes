#' ---------------------------
#'
#' Purpose of script:
#'
#' Compute IPSSM score using cross-validation with fold-specific training
#' Adapted from recomputeIPSSM.R to use cox_cv_functions.R
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' Input: results/gnn/preprocess/clinical_fold.tsv (with fold information)
#' Output: CV results by fold + aggregated results
#'
#' ---------------------------

## Load data and libraries
library(tidyverse)
library(ipssm)
library(survival)

# Load CV functions
source("scripts/3.GeneticScores/utils/cox_cv_functions.R")

## Load data - using clinical_fold.tsv as main input
clinical_train <- read_tsv("results/gnn/preprocess/clinical_fold.tsv")
# cyto <- read_table("data/IPSSMol/df_cna.tsv")
# mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
maf <- read.table("./data/IPSSMol/maf.tsv", header=T, sep="\t", stringsAsFactors = F)

# # Join data (clinical_fold.tsv should already have train info and folds)
# clinical_train <- clinical_raw %>%
#     left_join(cyto, by = "ID") %>%
#     left_join(mutation, by = "ID") 

## Modify variables to match IPSSM names (same preprocessing as original)
clinical_train$del7_7q <- as.numeric(clinical_train$del7q | clinical_train$del7)
clinical_train$del17_17p <- as.numeric(clinical_train$del17p | clinical_train$del17)
clinical_train$complex <- car::recode(clinical_train$complex, "'complex'='1'; 'non-complex'='0'")
clinical_train$CYTO_IPSSR <- car::recode(clinical_train$CYTO_IPSSR, "'Very-Good'='Very Good' ; 'Int'='Intermediate' ; 'Very-Poor'='Very Poor'")
clinical_train$TP53mut <- as.vector(sapply(clinical_train$ID, function(x) sum(maf$ID==x & maf$GENE=="TP53")))
clinical_train$TP53mut <- car::recode(clinical_train$TP53mut, " '2'='2 or more' ; '3'='2 or more' ")
clinical_train$TP53loh <- 0 
clinical_train$TP53loh[which(clinical_train$chr17 %in% c("cnloh","del"))] <- 1
clinical_train$TP53maxvaf <- NA
clinical_train$FLT3[which(clinical_train$FLT3_ITD==1)] <- 1

## Process data for model (same as original)
clin.process <- IPSSMprocess(clinical_train) %>%
    mutate(nRes2 = nRes2mean)

model_vars <- names(formals(IPSSMmain)$betaValues)
model_vars <- model_vars[model_vars != ""]
clinical_train_vars <- clin.process %>%
    select(all_of(model_vars), "AGE", "MDS_TYPE", "SEX", "LFS_YEARS", "LFS_STATUS", "IPSSM_SCORE", "ID", "fold")  %>%
    as_tibble() %>%
    filter(!is.na(LFS_YEARS) & !is.na(LFS_STATUS))

# Verify fold column exists
if (!"fold" %in% names(clinical_train_vars)) {
    stop("Column 'fold' not found in data. Check that clinical_fold.tsv contains fold information.")
}

cat("Data prepared for CV:\n")
cat("- Total patients:", nrow(clinical_train_vars), "\n")
cat("- Unique folds:", length(unique(clinical_train_vars$fold)), "\n")
cat("- Model variables:", length(model_vars), "\n")

## Execute Cross-Validation using cox_cv_functions
cat("\nExecuting cross-validation...\n")

# Create data for CV (exclude metadata columns)
cv_data <- clinical_train_vars %>%
    select(all_of(model_vars), "LFS_YEARS", "LFS_STATUS", "fold")

# Run CV
cv_results <- run_cox_cv(cv_data, formula = model_vars, 
                        time_var = "LFS_YEARS", 
                        status_var = "LFS_STATUS", 
                        fold_col = "fold")

cat("CV completed:\n")
cat("- Global C-index:", round(cv_results$global_c_index, 4), "\n")
cat("- Mean C-index:", round(cv_results$mean_c_index, 4), "±", round(cv_results$sd_c_index, 4), "\n")
cat("- Number of folds:", cv_results$n_folds, "\n")

## Create output directory for fold results
dir.create("results/gnn/IPSSM_cv_folds", recursive = TRUE, showWarnings = FALSE)

## Save individual fold results
fold_predictions_all <- data.frame()

for (fold_name in names(cv_results$fold_results)) {
    fold_result <- cv_results$fold_results[[fold_name]]
    fold_id <- fold_result$fold_id
    
    # Get patients in this fold
    fold_patients <- clinical_train_vars[clinical_train_vars$fold == fold_id, ]
    
    # Create fold prediction data
    fold_pred_data <- fold_patients %>%
        select(ID, LFS_YEARS, LFS_STATUS, fold, IPSSM_SCORE) %>%
        mutate(
            IPSSM_SCORE_cv = fold_result$predictions,
            c_index_fold = fold_result$c_index,
            n_test = fold_result$n_test,
            n_events = fold_result$n_events
        )
    
    # Save individual fold
    fold_filename <- paste0("results/gnn/IPSSM_cv_folds/fold_", fold_id, "_predictions.tsv")
    write.table(fold_pred_data, file = fold_filename, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Accumulate for global file
    fold_predictions_all <- rbind(fold_predictions_all, fold_pred_data)
    
    cat("Saved fold", fold_id, ":", nrow(fold_pred_data), "patients\n")
}

## Save CV summary
cv_summary <- data.frame(
    fold = sapply(cv_results$fold_results, function(x) x$fold_id),
    c_index = sapply(cv_results$fold_results, function(x) x$c_index),
    n_test = sapply(cv_results$fold_results, function(x) x$n_test),
    n_events = sapply(cv_results$fold_results, function(x) x$n_events)
) %>%
    rbind(data.frame(
        fold = "GLOBAL",
        c_index = cv_results$global_c_index,
        n_test = sum(sapply(cv_results$fold_results, function(x) x$n_test)),
        n_events = sum(sapply(cv_results$fold_results, function(x) x$n_events))
    ))

write.table(cv_summary, file = "results/gnn/IPSSM_cv_folds/cv_summary.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

## Traditional analysis for comparison (same as original)
cat("\nExecuting traditional analysis for comparison...\n")

cox_df <- clinical_train_vars %>% 
    select(-IPSSM_SCORE, -ID, -fold)

all_cox_model <- coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ., data = cox_df)

## Train only on training set (if train column exists)
if ("train" %in% names(clinical_train_vars)) {
    cox_df_train <- clinical_train_vars %>% 
        filter(is.na(train) | train) %>%
        select(-train, -IPSSM_SCORE, -ID, -fold)
    
    train_cox_model <- coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ., data = cox_df_train)
    
    clinical_train_vars_pred <- clinical_train_vars %>%
        mutate(AGE = 0, SEX = "F", MDS_TYPE = "primary") %>%
        mutate(across(all_of(model_vars), function(x) x - mean(x, na.rm = TRUE))) %>%
        mutate(IPSSM_SCORE_new = log(predict(all_cox_model, newdata = ., type = "risk")),
               IPSSM_SCORE_train = log(predict(train_cox_model, newdata = ., type = "risk")))
    
    # Add CV predictions to comparison
    clinical_train_vars_pred <- clinical_train_vars_pred %>%
        left_join(fold_predictions_all %>% select(ID, IPSSM_SCORE_cv), by = "ID")
    
    cat("Correlations with original IPSSM:\n")
    cat("- All data model:", round(cor(clinical_train_vars_pred$IPSSM_SCORE_new, 
                                     clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs"), 4), "\n")
    cat("- Train only model:", round(cor(clinical_train_vars_pred$IPSSM_SCORE_train, 
                                       clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs"), 4), "\n")
    cat("- CV model:", round(cor(clinical_train_vars_pred$IPSSM_SCORE_cv, 
                              clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs"), 4), "\n")
    
} else {
    # No train column available
    clinical_train_vars_pred <- clinical_train_vars %>%
        mutate(AGE = 0, SEX = "F") %>%
        mutate(across(all_of(model_vars), function(x) x - mean(x, na.rm = TRUE))) %>%
        mutate(IPSSM_SCORE_new = log(predict(all_cox_model, newdata = ., type = "risk")))
    
    # Add CV predictions
    clinical_train_vars_pred <- clinical_train_vars_pred %>%
        left_join(fold_predictions_all %>% select(ID, IPSSM_SCORE_cv), by = "ID")
    
    cat("Correlations with original IPSSM:\n")
    cat("- All data model:", round(cor(clinical_train_vars_pred$IPSSM_SCORE_new, 
                                     clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs"), 4), "\n")
    cat("- CV model:", round(cor(clinical_train_vars_pred$IPSSM_SCORE_cv, 
                              clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs"), 4), "\n")
}

## Generate plots (same as original but with CV comparison)
dir.create("figures/gnn", recursive = TRUE, showWarnings = FALSE)

# Original plot
png("figures/gnn/recomputeIPSSM_cv.png", width = 800, height = 600)
plot(clinical_train_vars_pred$IPSSM_SCORE_new, clinical_train_vars_pred$IPSSM_SCORE,
     main = "Recomputed IPSSM vs Original",
     xlab = "Recomputed IPSSM (All Data)", ylab = "Original IPSSM")
abline(0, 1, col = "red")
dev.off()

# CV comparison plot
png("figures/gnn/recomputeIPSSM_cv_comparison.png", width = 1200, height = 600)
par(mfrow = c(1, 2))

plot(clinical_train_vars_pred$IPSSM_SCORE_new, clinical_train_vars_pred$IPSSM_SCORE,
     main = "All Data Model", xlab = "Recomputed IPSSM", ylab = "Original IPSSM")
abline(0, 1, col = "red")

plot(clinical_train_vars_pred$IPSSM_SCORE_cv, clinical_train_vars_pred$IPSSM_SCORE,
     main = "Cross-Validation Model", xlab = "CV IPSSM", ylab = "Original IPSSM")
abline(0, 1, col = "red")

dev.off()

# Detailed comparison plot with ggplot
if ("train" %in% names(clinical_train_vars_pred)) {
    png("figures/gnn/recomputeIPSSM_cv_detailed.png", width = 1200, height = 800)
    p1 <- clinical_train_vars_pred %>%
        mutate(train = ifelse(is.na(train) | train, "Train", "Test")) %>%
        ggplot(aes(x = IPSSM_SCORE_cv, y = IPSSM_SCORE)) +
        geom_point() +
        geom_smooth(method = "lm") +
        labs(title = "Cross-Validation IPSSM",
             x = "CV IPSSM", y = "Original IPSSM") +
        facet_grid(~ train) +
        theme_bw()
    
    print(p1)
    dev.off()
}

## Save final results
write.table(clinical_train_vars_pred, 
            file = "results/gnn/patient_recomputedIPSSM_cv.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save all fold predictions
write.table(fold_predictions_all, 
            file = "results/gnn/IPSSM_cv_folds/all_fold_predictions.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Cross-validation completed successfully!\n")
cat("Files generated:\n")
cat("- Individual fold predictions: results/gnn/IPSSM_cv_folds/fold_*_predictions.tsv\n")
cat("- CV summary: results/gnn/cv_folds/cv_summary.tsv\n")
cat("- All predictions combined: results/gnn/IPSSM_cv_folds/all_fold_predictions.tsv\n")
cat("- Main results with CV: results/gnn/patient_recomputedIPSSM_cv.tsv\n")
cat("- Plots: figures/gnn/recomputeIPSSM_cv*.png\n")
cat("\nGlobal CV C-index:", round(cv_results$global_c_index, 4), "\n")
