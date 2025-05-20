#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide of MDS clusters using hierarchical clustering
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.5 R
#'
#' ---------------------------



# Load libraries and data
library(MASS)
library(cluster)
library(tidyverse)


load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select variables for clustering
#' - Demographic: age, sex
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
demo_vars <- c("AGE", "Sex")
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p",
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

# Compute distance for each type
#' Convert sex to dummy
#' Scale age
#' Convert PLT to log
#' Scale clin_vars
#' Select complete cases
clust_dataset <- mutate(clinical, Sex = ifelse(SEX == "F", 1, 0),
        PLT = log(PLT)) %>% 
    filter(consensus == "Low blasts") %>%
    select(all_of(c("ID", "Sex", "AGE", clin_vars, kar_events, mut_vars))) %>%
    filter(complete.cases(.))

#' Age distance	
age_dist <- dist(scale(clust_dataset$AGE), method = "euclidean")

#' Sex distance
sex_dist <- dist(clust_dataset$Sex, method = "manhattan")

#' Clinical distance
clin_mat <- clust_dataset[, clin_vars] %>%
    scale()
clin_dist <- dist(clin_mat, method = "euclidean")
clin_dist_norm <- clin_dist/max(clin_dist)

#' Karyotype distance
kar_mat <- clust_dataset[, kar_events]
kar_dist <- dist(kar_mat, method = "manhattan")
kar_dist_norm <- kar_dist/max(kar_dist)

#' Mutation distance
mut_mat <- clust_dataset[, mut_vars]
mut_dist <- dist(mut_mat, method = "manhattan")
mut_dist_norm <- mut_dist/max(mut_dist)

#' Compute overall distance
#' Normalize per maximum distance
comb_dist <- clin_dist_norm + kar_dist_norm + mut_dist_norm

png("figures/lowblast_hc_cluster/dist_distr.png")
par(mfrow = c(4, 1))
hist(clin_dist_norm, main = "Clinical distance", xlab = "Distance")
hist(kar_dist_norm, main = "Karyotype distance", xlab = "Distance")
hist(mut_dist_norm, main = "Mutation distance", xlab = "Distance")
hist(comb_dist, main = "Combined distance", xlab = "Distance")
dev.off()


png("figures/lowblast_hc_cluster/dist_cor.png")
par(mfrow = c(3, 1))
plot(clin_dist_norm[1:10000], comb_dist[1:10000], main = "Comb vs Clinical distance", xlab = "Clinical", ylab = "Combined")
plot(kar_dist_norm[1:10000], comb_dist[1:10000], main = "Comb vs Karyotype distance", xlab = "Karyotype", ylab = "Combined")
plot(mut_dist_norm[1:10000], comb_dist[1:10000], main = "Comb vs Mutation distance", xlab = "Mutation", ylab = "Combined")
dev.off()


png("figures/lowblast_hc_cluster/dist_cor2.png")
par(mfrow = c(3, 1))
plot(clin_dist_norm[1:10000], kar_dist_norm[1:10000], main = "Clinical vs Karyotype distance", xlab = "Clinical", ylab = "Karyotype")
plot(clin_dist_norm[1:10000], mut_dist_norm[1:10000], main = "Clinical vs Mutation distance", xlab = "Clinical", ylab = "Mutation")
plot(kar_dist_norm[1:10000], mut_dist_norm[1:10000], main = "Karyotype vs Mutation distance", xlab = "Karyotype", ylab = "Mutation")
dev.off()


cor(clin_dist_norm, comb_dist)
cor(kar_dist_norm, comb_dist)
cor(mut_dist_norm, comb_dist)

cor(clin_dist_norm, kar_dist_norm)
cor(clin_dist_norm, mut_dist_norm)
cor(kar_dist_norm, mut_dist_norm)


## Run hierarchical clustering
clust_hc <- hclust(comb_dist)

## Compute silhouettes for each subtype
silhouette_score_tree <- function(data_dist, clusth, k) {
  km <- cutree(clusth, k = k)
  ss <- silhouette(km, data_dist)
  mean(ss[, 3])
}

sil_scores <- sapply(2:40, silhouette_score_tree,
     data_dist = comb_dist, clusth = clust_hc)


png("figures/lowblast_hc_cluster/comb_hc_silhouette.png")
plot(2:40, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Select optimal number of clusters
k_opt <- which.max(sil_scores) + 1
k_opt

## Cut tree
clust_cut <- cutree(clust_hc, k = k_opt)
names(clust_cut) <- clust_dataset$ID

## Distance PCA
dist_pca <- isoMDS(comb_dist, k = 5)

png("figures/lowblast_hc_cluster/dist_pcoa.png")
plot(dist_pca$points, pch = 19, col = clust_cut)
dev.off()




clin_hc <- hclust(clin_dist)


sil_scores_clin <- sapply(2:40, silhouette_score_tree,
     data_dist = clin_dist, clusth = clin_hc)


png("figures/lowblast_hc_cluster/clin_hc_silhouette.png")
plot(2:40, sil_scores_clin, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

clin_clust <- cutree(clin_hc, k = 6)

clust_dataset$clin_hcclust <- clin_clust

clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clust_dataset, aes(y = .data[[x]], x = factor(clin_hcclust))) +
    geom_boxplot() +
    theme_bw() 
})
png("figures/lowblast_hc_cluster/clinclust_clinical_vars.png")
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()