#' ---------------------------
#'
#' Purpose of script:
#'
#' Prepare mutation data from deep learning models
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.4 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)

load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
maf <- read_tsv("./data/IPSSMol/maf.tsv")

## Select individuals with complete OS and clinical values
sel_vars <- c("AGE", "SEX", "BM_BLAST", "WBC", "ANC", 
    "MONOCYTES", "HB", "PLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "OS_YEARS", "OS_STATUS", "ID")

clin_comp <- clinical %>%
    select(sel_vars) %>%
    filter(complete.cases(.)) %>%
    mutate(across(c(AGE, BM_BLAST,  WBC, ANC, MONOCYTES, HB, PLT), scale),
        SEX = ifelse(SEX == "M", 1, 0),
        complex = ifelse(complex == "complex", 1, 0))

## Create mapping from individuals to mutations
mutation_map <- mutation %>% 
    filter(ID %in% clin_comp$ID) %>%
    pivot_longer(cols = -1, names_to = "Gene") %>%
    filter(value == 1 & !Gene %in% c("FLT3_ITD", "MLL_PTD", 
        "TP53mono", "TP53multi", "TET2bi", "TET2other")) %>%
    select(-value) 
    

mutation_maf <- maf %>%
    filter(ID %in% clin_comp$ID) %>%
    select(ID, GENE, VAF) 


write.table(clin_comp, file = "results/mutations/patients_features_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mutation_map, file = "results/mutations/mutation_map_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mutation_maf, file = "results/mutations/mutation_maf_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)


## Define another data.frame for testing the GNN model
clin_test <- clinical %>%
    select(sel_vars, any_of(unique(mutation_map$Gene))) %>%
    filter(complete.cases(.)) %>%
    mutate(across(c(AGE, BM_BLAST,  WBC, ANC, MONOCYTES, HB, PLT), scale),
        SEX = ifelse(SEX == "M", 1, 0),
        complex = ifelse(complex == "complex", 1, 0))
write.table(clin_test, file = "results/mutations/patients_gnn_test.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

comp_cases <- clinical %>%
    select(sel_vars) %>%
    filter(complete.cases(.))
scale_factors <- comp_cases %>%
    summarise(across(c(AGE, BM_BLAST,  WBC, ANC, MONOCYTES, HB, PLT), list(mean = mean, sd = sd)))

summary(clin_test$AGE - (comp_cases$AGE - scale_factors$AGE_mean)/scale_factors$AGE_sd)
summary(clin_test$BM_BLAST - (comp_cases$BM_BLAST - scale_factors$BM_BLAST_mean)/scale_factors$BM_BLAST_sd)
summary(clin_test$WBC - (comp_cases$WBC - scale_factors$WBC_mean)/scale_factors$WBC_sd)
summary(clin_test$ANC - (comp_cases$ANC - scale_factors$ANC_mean)/scale_factors$ANC_sd)
summary(clin_test$MONOCYTES - (comp_cases$MONOCYTES - scale_factors$MONOCYTES_mean)/scale_factors$MONOCYTES_sd)
summary(clin_test$HB - (comp_cases$HB - scale_factors$HB_mean)/scale_factors$HB_sd)
summary(clin_test$PLT - (comp_cases$PLT - scale_factors$PLT_mean)/scale_factors$PLT_sd)

colnames(scale_factors)[grep("BM", colnames(scale_factors))] <- c("BMBLAST_mean", "BMBLAST_sd") 

scale_factors_df <- scale_factors %>%
    pivot_longer(cols = everything(), names_to = c("Variable", ".value"), names_sep = "_") %>%
    mutate(Variable = case_when(
        Variable == "BMBLAST" ~ "BM_BLAST",
        TRUE ~ Variable
    ))
write.table(scale_factors_df, file = "results/mutations/scale_factors_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)