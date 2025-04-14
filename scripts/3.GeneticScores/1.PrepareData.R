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
library(caret)
library(survival)

load("results/preprocess/clinical_preproc.Rdata")
load("results/preprocess/ipssm_clinical_only.Rdata")
load("results/mutations/go_gene_map.Rdata")
load("results/tsne_perp10.Rdata")


mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
maf <- read_tsv("./data/IPSSMol/maf.tsv")

go_genes <- unique(tab_out$gene)

## Create table setting gene value to VAF
## Remove genes not present in GO terms
vaf_df <- maf %>% 
    select(ID, GENE, VAF) %>%
    filter(GENE %in% go_genes) %>%
    distinct() %>%
    pivot_wider(names_from = GENE, values_from = VAF, 
        values_fn = max, values_fill = 0)

## Add rows of samples without mutations
miss_ids <- clinical_all$ID[!clinical_all$ID %in% vaf_df$ID]

miss_df <- data.frame(ID = miss_ids, 
    matrix(0, nrow = length(miss_ids), ncol = ncol(vaf_df) -1 ))
colnames(miss_df) <- colnames(vaf_df)
vaf_df <- rbind(vaf_df, miss_df)

## Create data.frame indicating mutated/not-mutated
mut_df <- vaf_df
mut_df[, 2:ncol(mut_df)] <- lapply(mut_df[, 2:ncol(mut_df)], function(x){
    ifelse(x > 0, 1, 0)
})


## Substitute missing IPSSM by mean
ipssm_impute <- mutate(ipssm_clinonly, 
    IPSSMscore_imp = ifelse(is.na(IPSSMscore), IPSSMscore_mean, IPSSMscore))

## Combine all data into one data.frame
vaf_df_comb <- left_join(select(ipssm_impute, ID, IPSSMscore_imp), 
    vaf_df, by = "ID") %>%
    left_join(select(clinical_all, ID, IPSSM_SCORE, OS_YEARS, OS_STATUS), 
        by = "ID") %>%
        tibble()
mut_df_comb <- left_join(select(ipssm_impute, ID, IPSSMscore_imp),  
    mut_df, by = "ID") %>%
    left_join(select(clinical_all, ID, IPSSM_SCORE,  OS_YEARS, OS_STATUS), 
        by = "ID") %>%
        tibble()
save(vaf_df_comb, mut_df_comb, file = "results/mutations/mutation_table.Rdata")


## Remove missing survival data
miss_events <- is.na(mut_df_comb$OS_STATUS)

muts_comp <- subset(mut_df_comb, !miss_events)
vaf_comp <- subset(vaf_df_comb, !miss_events)

## Create test and train
set.seed(27)
train_indices <- createDataPartition(Surv(muts_comp$OS_YEARS, muts_comp$OS_STATUS),
     p = 0.75, list = FALSE)[,1]

write.table(muts_comp[train_indices, ], file = "results/mutations/boolean_mutations_train.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muts_comp[-train_indices, ], file = "results/mutations/boolean_mutations_test.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

write.table(vaf_comp[train_indices, ], file = "results/mutations/vaf_mutations_train.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(vaf_comp[-train_indices, ], file = "results/mutations/vaf_mutations_test.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)

## Create test and train using cluster
mut_df_clust <- mutate(mut_df_comb, cluster = factor(clust10_tree))
vaf_df_clust <- mutate(vaf_df_comb, cluster = factor(clust10_tree))


muts_comp_clust <- subset(mut_df_clust, !miss_events)
vaf_comp_clust <- subset(vaf_df_clust, !miss_events)

clust_vars_mut <- predict(dummyVars(~ cluster, data = muts_comp_clust), muts_comp_clust)
clust_vars_vaf <- predict(dummyVars(~ cluster, data = vaf_comp_clust), vaf_comp_clust)


write.table(cbind(muts_comp_clust[train_indices, 1], clust_vars_mut[train_indices, ], muts_comp_clust[train_indices, -c(1:2)]), 
    file = "results/mutations/boolean_cluster_mutations_train.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(muts_comp_clust[-train_indices, 1], clust_vars_mut[-train_indices, ], muts_comp_clust[-train_indices, -c(1:2)]), 
    file = "results/mutations/boolean_cluster_mutations_test.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)


write.table(cbind(vaf_comp_clust[train_indices, 1], clust_vars_vaf[train_indices, ], vaf_comp_clust[train_indices, -c(1:2)]), 
    file = "results/mutations/vaf_cluster_mutations_train.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(vaf_comp_clust[-train_indices, 1], clust_vars_vaf[-train_indices, ], vaf_comp_clust[-train_indices, -c(1:2)]), 
    file = "results/mutations/vaf_cluster_mutations_test.tsv", 
    sep = "\t", quote = FALSE, row.names = FALSE)


