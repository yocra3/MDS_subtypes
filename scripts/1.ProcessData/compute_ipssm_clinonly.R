#' ---------------------------
#'
#' Purpose of script:
#'
#' Compute IPSSM score assuming no molecular data is present
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.3 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)
library(ipssm)

load("results/preprocess/clinical_preproc.Rdata")

## Read ipssm example data
dd <- read.csv(system.file("extdata", "IPSSMexample.csv", package = "ipssm"), header = TRUE)

## Select variables present in model
model_df <- clinical_all %>%
    mutate(del7_7q = del7q, del17_17p = del17p) %>%
        select(ID, HB, PLT, BM_BLAST, del5q, del7_7q, complex,
    CYTO_IPSSR, del17_17p) %>%
    mutate(CYTO_IPSSR = gsub("-", " ", CYTO_IPSSR),
            CYTO_IPSSR = ifelse(CYTO_IPSSR == "Int", "Intermediate", CYTO_IPSSR))

## Add molecular variables as 0
extra <- colnames(dd)[!colnames(dd) %in% colnames(model_df)]
extra_df <- cbind(tibble(ID = model_df$ID), 
        matrix(0, nrow = nrow(model_df), ncol = length(extra))) %>%
        as_tibble()
colnames(extra_df) <- c("ID", extra)
model_df2 <- left_join(model_df, extra_df, by = "ID")

model_proc <- IPSSMprocess(model_df2)
model_res <- IPSSMmain(model_proc)
model_annot <- IPSSMannotate(model_res)

ipssm_clinonly <- select(model_annot, ID, starts_with("IPSSM"))
save(ipssm_clinonly, file = "results/preprocess/ipssm_clinical_only.Rdata")