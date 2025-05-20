#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide MDS low blast using PCA, CA and k-means
#' 
#' ---------------------------
#'
#' Notes:
#' Compute PCA for clinical variables and CA for karyotype and mutation
#' data. Use VAF for mutation data. Combine the components explaning >80% of variance in a new
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
library(corrplot)


load("results/preprocess/clinical_preproc.Rdata")
load("results/preprocess/mutation_patho_vaf_df.Rdata")


mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p", "complex",
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

# Select data for clustering
#' Select low blasts
#' Complete cases
#' Convert PLT to log
clust_dataset <- mutate(clinical, PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    select(all_of(c("ID", clin_vars, kar_events, mut_vars))) %>%
    filter(complete.cases(.)) %>%
    mutate(complex = ifelse(complex == "complex", 1, 0)) 

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

# Mutation PCA
mut_mat <- clust_dataset[, mut_vars] %>% as.matrix()
mut_mat  <- mut_mat [, colMeans(mut_mat ) > 0.05] ## Remove low frequency mutations & TET2
mut_mat <- mut_mat[, colnames(mut_mat) != "TET2"]
mut_mat <- (mut_mat - 1) * -1
mut_pca <- ca(mut_mat)

# Combine components
# Select components explaining ~80% of variation
comb_pc_mat <- cbind(clin_pc$x[, 1:4], kar_pca$rowcoord[, 1:4], mut_pca$rowcoord[, 1:7])
comb_pc <- prcomp(comb_pc_mat)

png("figures/blasts_filter_pc_cluster/corr_PCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat), method = 'number')
dev.off()

png("figures/blasts_filter_pc_cluster/corr_oriPCs_combPCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat, comb_pc$x[, 1:11]), method = 'number')
dev.off()


## K-means clustering
### Define number of clusters based on silhouette score
### Select 17 components (81.5% of variance)
set.seed(27)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

comb_pc_sel <- comb_pc$x[, 1:11]
dist_comb_pc <- dist(comb_pc_sel, method = "euclidean")
sil_scores <- sapply(2:20, silhouette_score,
      mat = comb_pc_sel, data_dist = dist_comb_pc)

png("figures/blasts_filter_pc_cluster/comb_kmeans_silhouette.png")
plot(2:20, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 7 and 16
clust_comb <- kmeans(comb_pc_sel, centers = 7, nstart = 1000)$cluster

tree_dataset <- select(clust_dataset, 
    all_of(c(colnames(clinvars_mat), colnames(kar_mat), colnames(mut_mat)))) %>%
    mutate(clust = clust_comb)

tree_mod <- rpart(clust ~ . , data = tree_dataset, method = "class")
predictions <- predict(tree_mod, tree_dataset, type = "class")
table(predictions, clust_comb)

png("figures/blasts_filter_pc_cluster/blasts_classification_tree_ori.png", width = 1000, height = 650)
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

png("figures/blasts_filter_pc_cluster/blasts_classification_tree.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()

## Relabel clusters and simplify classification
labels <- c("STAG2/8+", "del20q", "Y-", "RUNX1", "Highly Leukopenic", "Mildly Leukopenic", "7-" )
tree_data_slim <- tree_dataset %>% 
    mutate(clust = factor(labels[as.numeric(clust)], 
           levels = c("Highly Leukopenic", "Mildly Leukopenic", "Y-", "7-", "del20q", "RUNX1", "STAG2/8+"))) %>%
    select("WBC", "RUNX1", "STAG2", "del20q", "plus8", "del7", "delY", clust)

tree_labels <- rpart(clust ~ . , data = tree_data_slim, method = "class")


png("figures/blasts_filter_pc_cluster/blasts_classification_tree_labels.png", width = 1000, height = 650)
plot(tree_labels)
text(tree_labels, use.n = TRUE)
dev.off()
table(new = predict(tree_labels, tree_data_slim, type = "class"), ori = tree_data_slim$clust)


save(sil_scores, clin_pc, kar_pca, mut_pca, comb_pc, file = "results/clustering/blast_filter_pc_clustering_raw.Rdata")


clust_dataset$clust <-  predict(tree_labels, tree_dataset, type = "class")

load("results/preprocess/clinical_c_index.Rdata")

clinical_blasts <- left_join(clust_dataset, 
    select(clinical, hma, transplant, consensus, "ID", "IPSSM", "IPSSM_SCORE", "IPSSR", "IPSSR_SCORE", "AGE", "SEX", complex, starts_with("OS"), starts_with("AMLt")), by = "ID") %>%
    left_join(select(clinical_c_index, ID, OV_C_INDEX, SUB_C_INDEX), by ="ID")
save(tree_labels, clinical_blasts, file = "results/clustering/blast_filter_pc_clustering.Rdata")




## Define 16 clusters
clust_comb16 <- kmeans(comb_pc_sel, centers = 16, nstart = 1000)$cluster

lapply(c(colnames(kar_mat), colnames(mut_mat)), function(var){
    print(var)
    print(table(clust_dataset[[var]], clust_comb16))
})

tree_dataset16 <- select(clust_dataset, 
    all_of(c(colnames(clinvars_mat), colnames(kar_mat), colnames(mut_mat)))) %>%
    mutate(clust = clust_comb16)

tree_mod16 <- rpart(clust ~ . , data = tree_dataset16, method = "class")
predictions_16 <- predict(tree_mod16, tree_dataset16, type = "class")
table(predictions_16, clust_comb16)

png("figures/blasts_filter_pc_cluster/blasts_classification_tree_16c.png", width = 1000, height = 650)
plot(tree_mod16)
text(tree_mod16, use.n = TRUE)
dev.off()

clust11_df <- mutate(clust_dataset,
    clust11 = clust_comb11)
save(clust11_df, file = "results/clustering/lowblast_filter_pc_clustering_11C.Rdata")

## Re-cluster remaining samples
comb_pc_mat_mini <- comb_pc_mat[predictions %in% c(2, 3),]
comb_pc_mini <- prcomp(comb_pc_mat_mini, scale = TRUE)

comb_pc_mini_sel <- comb_pc_mini$x[, 1:13]
dist_comb_pc_mini <- dist(comb_pc_mini_sel, method = "euclidean")
set.seed(27)
sil_scores_mini <- sapply(2:20, silhouette_score,
      mat = comb_pc_mini_sel, data_dist = dist_comb_pc_mini)

png("figures/blasts_filter_pc_cluster/comb_kmeans_silhouette_mini.png")
plot(2:20, sil_scores_mini, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()


clust_comb <- kmeans(comb_pc_mat, centers = 5, nstart = 1000)$cluster


## Explore the effect of K on clustering
ks <- 2:11
clinical_low$complex <- factor(clinical_low$complex, levels = c("non-complex", "complex"))

lapply(ks, function(k){
    clust_comb <- kmeans(comb_pc$x[, 1:8], centers = k, nstart = 1000)$cluster
    print(table(clust_comb, clinical_low$clust))
    kars <- c(kar_events, "complex")
    kar_tabs <- lapply(kars, function(x) table(clinical_low[[x]], clust_comb))
    a <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
    colnames(a) <- kars
    print(a)

})
