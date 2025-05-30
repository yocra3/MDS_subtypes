#' ---------------------------
#'
#' Purpose of script:
#'
#'  Justify the division of MDS into 11 subtypes
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for subgroups definition.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(MASS)
library(tidyverse)
library(corrplot)
library(pheatmap)

load("results/clustering/combined_pc_clustering.Rdata")
load("results/clustering/combined_pc_clustering_raw.Rdata")

## Define color palette
colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")

scale <- scale_color_manual(values = colors, name = "Sub-group")

## Combined PCA plot
combined_full$clust_manual <- factor(combined_full$clust_manual, 
 levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", 
    "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))
df_comb$sub_group <- factor(combined_full$clust_manual, 
    levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", 
    "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))
png("figures/GESMD_IWS_clustering/subgroup_definition/combined_PCs_subgroups.png", width = 1800, height = 1500, res = 300)
df_comb %>%
    ggplot(aes(x = PC1, y = PC2, color = sub_group)) +
    geom_point() +
    scale +
    xlab("PC1 (11.7%)") +
    ylab("PC2 (7.5%)") +
    theme_bw()
dev.off()


## Features proportions
combined_full <- mutate(combined_full,
    complex = ifelse(complex.y %in% c("complex", 1), 1, 0))
features <- c("TET2other", "TET2bi", "delY",  "plus8",  "del7", "del20q", "del7q", "complex", "STAG2" )
names(features) <- features
feat_tabs <- lapply(features, function(x) table(combined_full[[x]], combined_full$clust_manual))
feat_tabs_props <- sapply(feat_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/GESMD_IWS_clustering/subgroup_definition/feat_proportion.png", width = 1800, height = 1500, res = 300)
pheatmap(mat = feat_tabs_props, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE)
dev.off()



## PCs exploration
comb_pc_mat <- cbind(clin_pc$x[, 1:4], kar_pca$rowcoord[, 1:5], mut_pca$rowcoord[, 1:8])
colnames(comb_pc_mat) <- c("clin_PC1", "clin_PC2", "clin_PC3", "clin_PC4",
    "kar_CA1", "kar_CA2", "kar_CA3", "kar_CA4", "kar_CA5", "mut_CA1", "mut_CA2",
    "mut_CA3", "mut_CA4", "mut_CA5", "mut_CA6", "mut_CA7", "mut_CA8")

png("figures/GESMD_IWS_clustering/subgroup_definition/corr_PCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat), method = 'number')
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_definition/corr_oriPCs_combPCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat, comb_pc$x[, 1:13]), method = 'number')
dev.off()


## Silhouette score
png("figures/GESMD_IWS_clustering/subgroup_definition/silhouette_score.png", width = 400, height = 300)
data.frame(N_clusters = 2:20, Score = sil_scores) %>%
  mutate(color = ifelse(N_clusters == 11, "rest", "selected")) %>%
  ggplot(aes(x = N_clusters, y = Score)) +
  geom_point(aes(color = color)) +
  geom_line() +
  scale_color_manual(values = c("red", "black")) +
  xlab("Number of clusters") +
  ylab("Silhouette Score") +
  scale_x_continuous(breaks = seq(2, 20, 4)) +
  theme_bw() +
  theme(legend.position = "none") 
dev.off()

## Match manual clusters with original clusters
tab <- table(combined_full$clust_manual, combined_full$clust_ori)
tab <- tab[c("Highly Leukopenic", "Midly Leukopenic", "del7q", "8+", "Y-",
    "complex", "TET2 monoallelic", "STAG2", "del20q", "7-", "TET2 bi-allelic"), ]
write.table(tab, file = "results/GESMD_IWS_clustering/cluster_mapping_table.txt", 
    sep = "\t", quote = FALSE, col.names = NA)