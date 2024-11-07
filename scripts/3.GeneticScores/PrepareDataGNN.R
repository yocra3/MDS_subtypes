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
sel_vars <- c("AGE", "SEX", "BM_BLAST", "PB_BLAST", "WBC", "ANC", 
    "MONOCYTES", "HB", "PLT", "del5q", "complex", "OS_YEARS", "OS_STATUS", "ID")

clin_comp <- clinical_all %>%
    select(sel_vars) %>%
    filter(complete.cases(.)) %>%
    mutate(across(c(AGE, BM_BLAST, PB_BLAST, WBC, ANC, MONOCYTES, HB, PLT), scale),
        SEX = ifelse(SEX == "M", 1, 0),
        complex = ifelse(complex == "complex", 1, 0))

## Create mapping from individuals to mutations
mutation_map <- mutation %>% 
    subset(ID %in% clin_comp$ID) %>%
    pivot_longer(cols = -1, names_to = "Gene") %>%
    filter(value == 1 & !Gene %in% c("FLT3_ITD", "MLL_PTD", 
        "TP53mono", "TP53multi", "TET2bi", "TET2other")) %>%
    select(-value) %>%
    

write.table(clin_comp, file = "results/mutations/patients_features_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mutation_map, file = "results/mutations/mutation_map_gnn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)