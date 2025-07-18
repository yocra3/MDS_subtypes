
library(tidyverse)
library(survival)

## Load data
clinical_train <- read_tsv("results/gnn/preprocess/clinical_fold.tsv")
recomputed_IPSSM <- read_table("results/gnn/patient_recomputedIPSSM_cv.tsv")
gnn_predictions <- read_tsv("results/gnn/full_input_model//combined_cv_predictions.tsv")

## Merge prediction with survival data
c_index_df <- clinical_train %>%
  select(ID, fold, OS_YEARS, OS_STATUS, AMLt_YEARS, AMLt_STATUS) %>%
  left_join(recomputed_IPSSM %>% select(ID, IPSSM_SCORE_cv), by = "ID") %>%
  left_join(gnn_predictions, by = "ID")

models <- colnames(c_index_df)[7:15]

compute_c_index_fold <- function(df, fold_id, model) {
  df_fold <- df %>% filter(fold == fold_id)
  
  # Compute C-index for OS
  surv_obj <- Surv(df_fold$OS_YEARS, df_fold$OS_STATUS)
  pred <- -df_fold[[model]]  # Invert for concordance
  c_index_os <- concordance(surv_obj ~ pred)$concordance
  
  # Compute C-index for AMLt
  surv_obj_amlt <- Surv(df_fold$AMLt_YEARS, df_fold$AMLt_STATUS)
  pred_amlt <- -df_fold[[model]]  # Invert for concordance
  c_index_amlt <- concordance(surv_obj_amlt ~ pred_amlt)$concordance

  return(data.frame(fold = fold_id, model = model, c_index_os = c_index_os, c_index_amlt = c_index_amlt))
}

c_index_sum <- lapply(unique(c_index_df$fold), function(fold) {
  lapply(models, function(model) {
    compute_c_index_fold(c_index_df, fold, model)
  }) %>%
    bind_rows() %>%
    mutate(fold = fold) 
}) %>%
  bind_rows() %>%
  arrange(fold)


c_index_sum %>% group_by(model) %>%
  summarise(mean_c_index_os = mean(c_index_os, na.rm = TRUE),
            sd_c_index_os = sd(c_index_os, na.rm = TRUE),
            mean_c_index_amlt = mean(c_index_amlt, na.rm = TRUE),
            sd_c_index_amlt = sd(c_index_amlt, na.rm = TRUE)) %>%
  write_tsv("results/gnn/full_input_model/c_index_OS_AMLt_summary.tsv")

