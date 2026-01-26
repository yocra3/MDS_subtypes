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
library(FactoMineR)

load("results/clustering/MFA_results.Rdata")


## Define color palette
# colors11 <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
#     "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")
colors6 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#999999", "grey40",  "black")
scale <- scale_color_manual(values = colors6, name = "Sub-group")

classifySamples <- function(df){
    new_class <- case_when(
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    df$del7 == 1       ~ "7-",
    df$STAG2 == 1       ~ "STAG2",
    df$BM_BLAST <= 5    ~ "Low blasts",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c( "EZH2", "TET2 bi-allelic", "7-", "STAG2", 
      "Low blasts", "MDS-IB1", "MDS-IB2"))
}

## MFA plots
###############################################################################

## IWS
mfa_df <- tibble(Sample = IWS_dataset$ID,
    PC1 = IWS_mfa$ind$coord[, 1],
    PC2 = IWS_mfa$ind$coord[, 2]) %>%
    left_join(IWS_dataset, by = c("Sample" = "ID")) %>%
    mutate(sub_group = classifySamples(.),
    oriCluster = factor(clust_comb1))

iws_mfa_plot <- mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = sub_group)) +
    geom_point() +
    ggtitle("IWS") +
    scale +
    xlab("Dim 1 (7.4%)") +
    ylab("Dim 2 (5.6%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA.png", width = 1500, height = 1200, res = 300)
iws_mfa_plot
dev.off()


png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_oriCluster.png", width = 1500, height = 1200, res = 300)
 mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = oriCluster)) +
    geom_point() +
    ggtitle("IWS") +
    xlab("Dim 1 (7.4%)") +
    ylab("Dim 2 (5.6%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()


## GESMD
gesmd_mfa_df <- tibble(Sample = gesmd_dataset_filt$ID,
    PC1 = GESMD_mfa$ind$coord[, 1],
    PC2 = GESMD_mfa$ind$coord[, 2]) %>%
    left_join(gesmd_dataset_filt, by = c("Sample" = "ID")) %>%
    mutate(sub_group = classifySamples(.),
    oriCluster = factor(clust_gesmd))

gesmd_mfa_plot <- gesmd_mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = sub_group)) +
    geom_point() +
    ggtitle("GESMD") +
    scale +
    xlab("Dim 1 (6.2%)") +
    ylab("Dim 2 (5.6%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA.png", width = 1500, height = 1200, res = 300)
gesmd_mfa_plot
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA_oriCluster.png", width = 1500, height = 1200, res = 300)
 gesmd_mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = oriCluster)) +
    geom_point() +
    ggtitle("GESMD") +
    xlab("Dim 1 (6.2%)") +
    ylab("Dim 2 (5.6%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Features proportions cluster
features_clustIWS <- c("TET2bi", "EZH2",  "PHF6",  "STAG2", "Other Functions")
names(features_clustIWS) <- features_clustIWS
feat_tabs_clustIWS <- lapply(features_clustIWS, function(x) table(IWS_cluster[[x]], clust_comb1))
feat_tabs_clustIWS[["7- or SETPBP1"]] <- table(IWS_cluster[["del7"]] == 1 | IWS_cluster[["SETBP1"]] == 1, clust_comb1)
feat_tabs_props_clustIWS <- sapply(feat_tabs_clustIWS, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_oricluster_feat_proportion.png", width = 1000, height = 1000, res = 300)
pheatmap(mat = feat_tabs_props_clustIWS, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "IWS", legend = FALSE)
dev.off()


features_clustGESMD <- c("TET2bi", "EZH2", "STAG2", "CBL", "RAS MAPK Pathway", "del7")
names(features_clustGESMD) <- features_clustGESMD
feat_tabs_clustGESMD <- lapply(features_clustGESMD, function(x) table(GESMD_cluster[[x]], clust_gesmd))
feat_tabs_clustGESMD[["U2AF1 or del20q"]] <- table(gesmd_mfa_df[["U2AF1"]] == 1 | gesmd_mfa_df[["del20q"]] == 1, clust_gesmd)
feat_tabs_props_clustGESMD <- sapply(feat_tabs_clustGESMD, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_oricluster_feat_proportion.png", width = 1200, height = 1300, res = 300)
pheatmap(mat = feat_tabs_props_clustGESMD, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "GESMD", legend = FALSE)
dev.off()
## Features proportions
features <- c("TET2bi", "EZH2",  "del7", "STAG2")
names(features) <- features
feat_tabs <- lapply(features, function(x) table(mfa_df[[x]], mfa_df$sub_group))
feat_tabs_props <- sapply(feat_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_feat_proportion.png", width = 1000, height = 1000, res = 300)
pheatmap(mat = feat_tabs_props, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "IWS", legend = FALSE)
dev.off()

gesmd_feat_tabs <- lapply(features, function(x) table(gesmd_mfa_df[[x]], gesmd_mfa_df$sub_group))
gesmd_feat_tabs_props <- sapply(gesmd_feat_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_feat_proportion.png", width = 1000, height = 1000, res = 300)
pheatmap(mat = gesmd_feat_tabs_props, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "GESMD", legend = FALSE)
dev.off()



## PCs exploration
IWS_dim_cors <- rbind(tibble(Dim1 = sqrt(IWS_mfa$quali.var$cos2[, 1]),
    Dim2 = sqrt(IWS_mfa$quali.var$cos2[, 2]),
    pos1 = IWS_mfa$quali.var$coord[, 1],
    pos2 = IWS_mfa$quali.var$coord[, 2], 
    Variable = rownames(IWS_mfa$quali.var$cos2)) %>%
    mutate(Dim1 = ifelse(pos1 < 0, -Dim1, Dim1),
        Dim2 = ifelse(pos2 < 0, -Dim2, Dim2)) %>%
    filter(grepl("_1", Variable)) %>%
    mutate(Variable = gsub("_1", "", Variable)) %>%
    filter(Variable %in% c("TET2bi", "EZH2",  "PHF6",  "del7", "STAG2")) %>%
    select(Dim1, Dim2, Variable),
    tibble(Dim1 = IWS_mfa$quanti.var$cor[, 1],
    Dim2 = IWS_mfa$quanti.var$cor[, 2],
    Variable = rownames(IWS_mfa$quanti.var$cor))
)

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_dimension_cors.png", width = 1800, height = 1500, res = 300)
ggplot(IWS_dim_cors, aes(x = Dim1, y = Dim2, label = Variable)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
        arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(vjust = -0.5) +
    xlab("Dim 1 correlation") +
    ylab("Dim 2 correlation") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_bw()
dev.off()

## Selected components
png("figures/GESMD_IWS_clustering/subgroup_definition/joint_components_inertia.png", width = 400, height = 300)
tibble(Dimension = seq_len(40), IWS = IWS_mfa$global.pca$eig[1:40, 2], GESMD = GESMD_mfa$global.pca$eig[1:40, 2]) %>%
    pivot_longer(cols = c("IWS", "GESMD"), names_to = "Dataset", values_to = "Inertia") %>%
    mutate(Selection = ifelse((Dimension <= 7 & Dataset == "IWS") | (Dimension <= 9 & Dataset == "GESMD"), "chosen", "rest"),
    Dataset = factor(Dataset, levels = c("IWS","GESMD")))  %>%
    ggplot(aes(x = Dimension, y = Inertia, fill = Selection)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("black", "white")) +
    facet_grid(. ~ Dataset) +
    xlab("Dimension") +
    ylab("Percentage of inertia (%)") +
    theme_bw() +
    theme(legend.position = "none")
dev.off()



## Silhouette score
png("figures/GESMD_IWS_clustering/subgroup_definition/joint_silhouette_score.png", width = 400, height = 300)
tibble(N_clusters = 2:20, IWS = sil_scores1, GESMD = sil_scores_gesmd) %>%
    pivot_longer(cols = c("IWS", "GESMD"), names_to = "Dataset", values_to = "Score") %>%
    mutate(color = ifelse((N_clusters == 8 & Dataset == "IWS") | (N_clusters == 10 & Dataset == "GESMD"), "selected", "rest")) %>%
  ggplot(aes(x = N_clusters, y = Score, color = Dataset)) +
  geom_point(aes(size = color)) +
  geom_line() +
  xlab("Number of clusters") +
  ylab("Silhouette Score") +
  scale_x_continuous(breaks = seq(2, 20, 4)) +
  theme_bw() +
  guides(size = "none") 
dev.off()

## Match manual clusters with original clusters
IWS_tab <- table(mfa_df$sub_group, clust_comb1)
IWS_tab <- IWS_tab[levels(mfa_df$sub_group), ]

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_cluster_mapping.png", width = 1800, height = 1500, res = 300)
pheatmap(mat = IWS_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0)
dev.off()

GESMD_tab <- table(gesmd_mfa_df$sub_group, clust_gesmd)
GESMD_tab <- GESMD_tab[levels(gesmd_mfa_df$sub_group), ]

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_cluster_mapping.png", width = 1800, height = 1500, res = 300)
pheatmap(mat = GESMD_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0)
dev.off()

