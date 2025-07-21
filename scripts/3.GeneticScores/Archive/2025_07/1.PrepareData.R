#' ---------------------------
#'
#' Purpose of script:
#'
#' Prepare data from deep learning models
#' We will generate two data.frames:
#' 1. One with clinical variables and karyotipic variables. This data.frame 
#' will feed the patient variables for the GNN model. Additionally,
#' we will include the mutations as boolean variables for training the NN models,
#' and the survival data and IPSSM scores to validate the model.
#' 2. One with the mutations, including the VAF. This data.frame will be used
#' to define the gene nodes in the GNN model.
#'
#' ---------------------------
#'
#' Notes:
#'  Do not include age into the model. Age is related to survival but may not be
#' related to disease. Exclude MONOCYTES due to the high number of missing values.
#' 
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.5 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)
library(caret)
library(survival)

load("results/preprocess/clinical_preproc.Rdata")

## Load mutation data and select relevant columns
maf <- read_tsv("./data/IPSSMol/maf.tsv") %>%
    select(ID, GENE, VAF)

## Select individuals with complete OS, clinical values, karyotype and mutations
model_vars <-  c("SEX", "BM_BLAST", "WBC", "ANC", 
    "HB", "PLT", "logPLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex",  "TP53multi",
    "NPM1", "RUNX1", "NRAS", "ETV6", "IDH2", "CBL", "EZH2", "U2AF1",
    "SRSF2", "DNMT3A", "ASXL1", "KRAS",  "ID", "LFS_YEARS", "LFS_STATUS")

sel_vars <- c(model_vars, "OS_YEARS", "OS_STATUS", "AMLt_YEARS", "AMLt_STATUS", "IPSSM_SCORE")

clinical$logPLT <- log(clinical$PLT)

comp_cases <- clinical %>%
    select(all_of(model_vars)) %>%
    filter(complete.cases(.))

clin_models <- clinical %>%
    select(all_of(sel_vars)) %>%
    filter(ID %in% comp_cases$ID) %>%
    mutate(across(c(BM_BLAST, WBC, ANC, HB, PLT, logPLT), scale),
        SEX = ifelse(SEX == "M", 1, 0), 
        complex = ifelse(complex == "complex", 1, 0))

## Create test and train
set.seed(27)
train_indices <- createDataPartition(Surv(clin_models$LFS_YEARS, clin_models$LFS_STATUS),
     p = 0.8, list = FALSE)[,1]
clin_models$train <- seq_len(nrow(clin_models)) %in% train_indices | is.na(clin_models$IPSSM_SCORE)

write.table(clin_models %>% select(-logPLT), file = "results/gnn/preprocess/patient_variables.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

write.table(clin_models %>% select(-PLT), file = "results/gnn/preprocess/patient_variables_logPLT.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)


write.table(maf, file = "results/gnn/preprocess/mutation_maf_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

## Recover scale factors
scale_factors <- comp_cases %>%
    summarise(across(c(BM_BLAST,  WBC, ANC, HB, PLT, logPLT), list(mean = mean, sd = sd)))

summary(clin_models$BM_BLAST - (comp_cases$BM_BLAST - scale_factors$BM_BLAST_mean)/scale_factors$BM_BLAST_sd)
summary(clin_models$WBC - (comp_cases$WBC - scale_factors$WBC_mean)/scale_factors$WBC_sd)
summary(clin_models$ANC - (comp_cases$ANC - scale_factors$ANC_mean)/scale_factors$ANC_sd)
summary(clin_models$HB - (comp_cases$HB - scale_factors$HB_mean)/scale_factors$HB_sd)
summary(clin_models$PLT - (comp_cases$PLT - scale_factors$PLT_mean)/scale_factors$PLT_sd)
summary(clin_models$logPLT - (comp_cases$logPLT - scale_factors$logPLT_mean)/scale_factors$logPLT_sd)
colnames(scale_factors)[grep("BM", colnames(scale_factors))] <- c("BMBLAST_mean", "BMBLAST_sd") 

scale_factors_df <- scale_factors %>%
    pivot_longer(cols = everything(), names_to = c("Variable", ".value"), names_sep = "_") %>%
    mutate(Variable = case_when(
        Variable == "BMBLAST" ~ "BM_BLAST",
        TRUE ~ Variable
    ))
write.table(scale_factors_df, file = "results/gnn/preprocess/scale_factors_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)