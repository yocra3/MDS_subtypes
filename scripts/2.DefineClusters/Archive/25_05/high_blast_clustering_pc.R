#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide MDS high blast using PCA, CA and k-means
#' 
#' ---------------------------
#'
#' Notes:
#' Compute PCA for clinical variables and CA for karyotype and mutation
#' data. Combine the components explaning >80% of variance in a new
#' matrix and run PCA again. Run k-means clustering on the final matrix.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(MASS)
library(cluster)
library(tidyverse)
library(rpart)
library(ca)


load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p",
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

# Select data for clustering
#' Select low blasts
#' Complete cases
#' Convert PLT to log
clust_dataset <- mutate(clinical, PLT = log(PLT)) %>% 
    filter(consensus %in% c("MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
    select(all_of(c("ID", clin_vars, kar_events, mut_vars))) %>%
    filter(complete.cases(.))

# Clinical PCA
clinvars_mat <- clust_dataset %>%
    select(BM_BLAST, PB_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
    as.matrix()

clin_pc <- prcomp(clinvars_mat, scale = TRUE)

# Karyotype CA
kar_mat <- clust_dataset[, kar_events] %>% as.matrix()
kar_mat  <- kar_mat [, colMeans(kar_mat ) > 0.01] ## Remove low frequency karyotypes
kar_mat <- (kar_mat - 1) * -1
kar_pca <- ca(kar_mat)

# Mutation CA
mut_mat <- clust_dataset[, mut_vars] %>% as.matrix()
mut_mat  <- mut_mat [, colMeans(mut_mat ) > 0.05] ## Remove low frequency mutations
mut_mat <- (mut_mat - 1) * -1
mut_pca <- ca(mut_mat)

# Combine components
# Select components explaining ~80% of variation
comb_pc_mat <- cbind(clin_pc$x[, 1:4], kar_pca$rowcoord[, 1:4], mut_pca$rowcoord[, 1:10])
comb_pc <- prcomp(comb_pc_mat)

## K-means clustering
### Define number of clusters based on silhouette score
### Select 8 components (76.8% of variance)
set.seed(27)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

comb_pc_sel <- comb_pc$x[, 1:13]
dist_comb_pc <- dist(comb_pc_sel, method = "euclidean")
sil_scores <- sapply(2:20, silhouette_score,
      mat = comb_pc_sel, data_dist = dist_comb_pc)



png("figures/high_blasts_filter_pc_cluster/comb_kmeans_silhouette.png")
plot(2:20, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()
clust_comb <- kmeans(comb_pc$x[, 1:13], centers = 7, nstart = 1000)$cluster

tree_dataset <- select(clust_dataset, 
    all_of(c(colnames(clinvars_mat), colnames(kar_mat), colnames(mut_mat)))) %>%
    mutate(clust = clust_comb)

tree_mod <- rpart(clust ~ . , data = tree_dataset, method = "class")
predictions <- predict(tree_mod, tree_dataset, type = "class")
table(predictions, clust_comb)

png("figures/high_blasts_filter_pc_cluster/high_blasts_classification_tree_ori.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()

## Re-assing samples to clusters
round <- 1
tree_class <- clust_comb
while (sum(tree_class != predictions) > 0) {
    tree_class <- predictions
    tree_dataset$clust <- predictions
    tree_mod <- rpart(clust ~ . , data = tree_dataset, method = "class")
    predictions <- predict(tree_mod, tree_dataset, type = "class")
    round <- round + 1
    print(round)
}

png("figures/high_blasts_filter_pc_cluster/high_blasts_classification_tree.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()

## Relabel clusters and simplify classification
labels <- c("7-", "High Leukopenic", "del20q", "complex", "Mild Leukopenic", "TET2 mono", "TET2 bi" )
tree_data_slim <- tree_dataset %>% 
    mutate(clust = factor(labels[as.numeric(clust)], 
           levels = c("High Leukopenic", "Mild Leukopenic", "TET2 mono", "TET2 bi", "del20q", "7-", "complex"))) %>%
    select( "WBC", "del7", "del5q", "TP53", "del20q", "TET2bi", "TET2", "ANC", clust)

tree_labels <- rpart(clust ~ . , data = tree_data_slim, method = "class")

png("figures/high_blasts_filter_pc_cluster/high_blasts_classification_tree_labels.png", width = 1000, height = 650)
plot(tree_labels)
text(tree_labels, use.n = TRUE)
dev.off()
table(new = predict(tree_labels, tree_data_slim, type = "class"), ori = tree_data_slim$clust)


save(sil_scores, tree_labels, clin_pc, kar_pca, mut_pca, comb_pc, file = "results/clustering/high_blast_filter_pc_clustering_raw.Rdata")


clust_dataset$clust <-  predict(tree_labels, tree_dataset, type = "class")

load("results/preprocess/clinical_c_index.Rdata")

clinical_blasts <- left_join(clust_dataset, 
    select(clinical, hma, transplant, consensus, "ID", "IPSSM", "IPSSM_SCORE", "IPSSR", "IPSSR_SCORE", "AGE", "SEX", complex, starts_with("OS"), starts_with("AMLt")), by = "ID") %>%
    left_join(select(clinical_c_index, ID, OV_C_INDEX, SUB_C_INDEX), by ="ID")
save(clinical_blasts, file = "results/clustering/high_blast_filter_pc_clustering.Rdata")

sil <- silhouette(clust_comb, dist = dist_comb_pc)
plot(sil, main = "hihg blast clusters",
     do.col.sort=FALSE)