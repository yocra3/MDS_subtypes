#' ---------------------------
#'
#' Purpose of script:
#'
#'  Get mutation data for passenger mutations
#' 
#' ---------------------------
#'
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(tidyverse)

load("results/preprocess/clinical_preproc.Rdata")
all_mutations <- read_delim("data/cbioportal/mds_iwg_2022/data_mutations.txt", delim = "\t") %>%
    mutate(ID = Tumor_Sample_Barcode,
    VAF = t_alt_count/(t_alt_count + t_ref_count ))
mutations_clin <- subset(all_mutations, ID %in% clinical$ID)

## Create VAF table with all mutations

## Add rows of samples without mutations

vaf_all_df <- mutations_clin %>% 
    select(ID, Hugo_Symbol, VAF) %>%
    distinct() %>%
    pivot_wider(names_from = Hugo_Symbol, values_from = VAF, 
        values_fn = max, values_fill = 0)

miss_ids <- clinical$ID[!clinical$ID %in% mutations_clin$ID]

miss_df <- data.frame(ID = miss_ids, 
    matrix(0, nrow = length(miss_ids), ncol = ncol(vaf_all_df) -1 ))
colnames(miss_df) <- colnames(vaf_all_df)
vaf_all_df <- rbind(vaf_all_df, miss_df)

save(vaf_all_df, file = "results/preprocess/mutation_all_vaf_df.Rdata")

## Create VAF table with pathogenic mutations
vaf_patho_df <- mutations_clin %>% 
    filter(cbp_driver == "Putative_Driver" ) %>%
    select(ID, Hugo_Symbol, VAF) %>%
    distinct() %>%
    pivot_wider(names_from = Hugo_Symbol, values_from = VAF, 
        values_fn = max, values_fill = 0)

miss_ids2 <- clinical$ID[!clinical$ID %in% vaf_patho_df$ID]

miss_df2 <- data.frame(ID = miss_ids2, 
    matrix(0, nrow = length(miss_ids2), ncol = ncol(vaf_patho_df) -1 ))
colnames(miss_df2) <- colnames(vaf_patho_df)
vaf_patho_df <- rbind(vaf_patho_df, miss_df2)

save(vaf_patho_df, file = "results/preprocess/mutation_patho_vaf_df.Rdata")


# Mutation PCA (VAF)
vaf_mat <- vaf_all_df[, -1] %>% as.matrix()
rownames(vaf_mat) <- vaf_all_df$ID
vaf_pca <- prcomp(vaf_mat)

vaf_mat2 <- vaf_patho_df[, -1] %>% as.matrix()
rownames(vaf_mat2) <- vaf_patho_df$ID
vaf_mat2 <- vaf_mat2[rownames(vaf_mat), ]
vaf_pca2 <- prcomp(vaf_mat2)
