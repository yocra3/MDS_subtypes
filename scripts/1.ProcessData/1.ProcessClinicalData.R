#' ---------------------------
#'
#' Purpose of script:
#'
#'  Process clinical data used for training IPSS molecular (v2 after talking with Irene and Teresa)
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.4 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
clinical_raw <- read_tsv("./data/IPSSMol/df_clinical.tsv")
cyto <- read_table("data/IPSSMol/df_cna.tsv")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

#' ---------------------------
#' Cytogenetics
# Remove missing events before merging
ev_freq <- colSums(cyto[, -1])
cyto_noNull <- cyto[, c(TRUE, ev_freq > 0)]

# Edit variables
#' Include columns to indicate if the number of aberrations, 
#' deletion, duplications, monosomy or trisomy
cyto_sum <- cyto_noNull %>%
  mutate(N_aberrations = rowSums(. == 1), 
         N_deletions = rowSums(select(., matches("del")) == 1),
         N_duplications = rowSums(select(., matches("plus")) == 1),
         N_rearrangements =  rowSums(select(., matches("r_")) == 1),
         N_monosomies = rowSums(select(., matches("del[0-9]$")) == 1),
         N_trisomies = rowSums(select(., matches("plus[0-9]$")) == 1))


#' Mutations
## Select mutations present in at least 10 samples
sel_muts <- colSums(mutation[, -1]) > 10 
sel_muts <- colnames(mutation[, -1])[sel_muts]

mutation_N <- mutate(mutation, 
    N_mutations = rowSums(mutation[, 2:127])) %>%
    select(sel_muts, ID, N_mutations)

clin_comb <- left_join(clinical_raw, cyto_sum, by = "ID") %>%
  left_join(mutation_N, by = "ID")
 
# Data transform
#' Create levels for different variables
#' Round clinical variables 
#' Consider that samples with missing RINGED_SIDEROBLASTS have 0 sideroblasts.
clinical <- mutate(clin_comb, 
    CYTO_IPSSR = factor(CYTO_IPSSR, levels = c("Very-Poor", "Poor", "Int", "Good", "Very-Good")),
    IPSSR = factor(IPSSR, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSRA = factor(IPSSRA, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
    OS_YEARS = ifelse(OS_YEARS == 0, 0.0001, OS_YEARS), ## Remove 0s in survival
    BM_BLAST = round(BM_BLAST, 0),
    PB_BLAST = round(PB_BLAST, 0),
    WBC = round(WBC, 0),
    ANC = round(ANC, 0), 
    MONOCYTES = round(MONOCYTES, 0),
    HB = round(HB, 0), 
    PLT = round(PLT, 0),
    RINGED_SIDEROBLASTS = ifelse(is.na(RINGED_SIDEROBLASTS), 0, RINGED_SIDEROBLASTS), 
    ## Compute consensus MDS sub-types
    consensus = ifelse(TP53multi == 1 & BM_BLAST <= 20, "Mutated TP53",
      ifelse(del5q == 1 & del7q == 0 & BM_BLAST <= 5, "del5q",
        ifelse(SF3B1 > 0 & del7q == 0 & complex == "non-complex" & BM_BLAST <= 5, "mutated SF3B1",
          ifelse(BM_BLAST <= 5, "Low blasts",
            ifelse(BM_BLAST > 10, "MDS-IB2",
              ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other"))))))
) %>%
  filter(!is.na(consensus)) ## Select samples with consensus classification
 
#' Merge all data.frames

save(clinical, file = "results/preprocess/clinical_preproc.Rdata")

