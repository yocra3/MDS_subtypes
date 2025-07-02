#' ---------------------------
#'
#' Purpose of script:
#'
#' Prepare data for NN models
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)

load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select individuals with complete OS, clinical values, karyotype and mutations
sel_vars <- c("AGE", "SEX", "BM_BLAST", "WBC", "ANC", 
    "MONOCYTES", "HB", "PLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "OS_YEARS", "OS_STATUS", "ID", "TP53multi",
    "NPM1", "RUNX1", "NRAS", "ETV6", "IDH2", "CBL", "EZH2", "U2AF1",
    "SRSF2", "DNMT3A", "ASXL1", "KRAS")

clin_NN <- clinical %>%
    select(all_of(sel_vars)) %>%
    filter(complete.cases(.)) %>%
    mutate(across(c(AGE, BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT), scale),
        SEX = ifelse(SEX == "M", 1, 0),


write.table(clin_NN, file = "results/mutations/patients_features_nn.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
