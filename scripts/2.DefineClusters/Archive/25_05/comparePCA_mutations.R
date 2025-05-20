#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compare PCA and CA for mutation data
#' 
#' ---------------------------
#'
#' Notes:
#' Compute running a PCA con mutation data based on VAF with a CA in
#' mutation data based on presence/absence of mutations.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)
library(ca)


load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
maf <- read_tsv("./data/IPSSMol/maf.tsv")

## Create table setting gene value to VAF
## Remove genes not present in GO terms
vaf_df <- maf %>% 
    select(ID, GENE, VAF) %>%
    distinct() %>%
    pivot_wider(names_from = GENE, values_from = VAF, 
        values_fn = max, values_fill = 0)

## Add rows of samples without mutations
miss_ids <- mutation$ID[!mutation$ID %in% vaf_df$ID]

miss_df <- data.frame(ID = miss_ids, 
    matrix(0, nrow = length(miss_ids), ncol = ncol(vaf_df) -1 ))
colnames(miss_df) <- colnames(vaf_df)
vaf_df <- rbind(vaf_df, miss_df)

save(vaf_df, file = "results/preprocess/mutation_vaf_df.Rdata")

# Mutation CA
mut_mat <- mutation[, -1] %>% as.matrix()
mut_mat <- (mut_mat - 1) * -1
rownames(mut_mat) <- mutation$ID
mut_pca <- ca(mut_mat)

# Mutation PCA (VAF)
vaf_mat <- vaf_df[, -1] %>% as.matrix()
rownames(vaf_mat) <- vaf_df$ID
vaf_mat <- vaf_mat[rownames(mut_mat), ]
vaf_mat[is.na(vaf_mat)] <- 0
vaf_pca <- prcomp(vaf_mat)
