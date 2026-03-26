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
library(cowplot)

load("results/clustering/MFA_results.Rdata")


## Define color palette
# colors11 <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
#     "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
    "#D55E00", "#999999", "grey40",  "black")
scale <- scale_color_manual(values = colors, name = "Sub-group")

classifySamples <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "Complex",
    df$del5q == 1      ~ "del5q-IB",
    df$SF3B1 == 1      ~ "SF3B1-IB",
    df$del7 == 1       ~ "-7",
    df$EZH2 == 1        ~ "EZH2",
    df$STAG2 == 1       ~ "STAG2",
    df$TET2bi == 1      ~ "TET2-bi",
    df$BM_BLAST <= 5    ~ "MDS-LB",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c("EZH2", "TET2-bi",  "-7", "STAG2", "del5q-IB", "SF3B1-IB",
   "Complex", "MDS-LB", "MDS-IB1", "MDS-IB2"))
}


## MFA plots
###############################################################################

## IWS
mfa_df <- tibble(Sample = IWS_dataset$ID,
    PC1 = IWS_mfa$ind$coord[, 1],
    PC2 = IWS_mfa$ind$coord[, 2],
    PC3 = IWS_mfa$ind$coord[, 3],
    PC4 = IWS_mfa$ind$coord[, 4]) %>%
    left_join(IWS_dataset, by = c("Sample" = "ID")) %>%
    mutate(sub_group = classifySamples(.),
    oriCluster = factor(clust_comb1))

iws_mfa_plot <- mfa_df %>%
    mutate(alpha = ifelse(sub_group %in% c("MDS-LB", "MDS-IB1", "MDS-IB2"), 0.3, 1)) %>%
    ggplot(aes(x = PC1, y = PC2, color = sub_group, alpha = alpha)) +
    geom_point() +
    scale_alpha_identity() +
    # ggtitle("IWS") +
    scale +
    xlab("Dim 1 (7.2%)") +
    ylab("Dim 2 (5.5%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA.png", width = 1500, height = 1200, res = 300)
iws_mfa_plot
dev.off()


iws_mfa_plot_d34 <- mfa_df %>%
    mutate(alpha = ifelse(sub_group %in% c("MDS-LB", "MDS-IB1", "MDS-IB2"), 0.3, 1)) %>%
    ggplot(aes(x = PC3, y = PC4, color = sub_group, alpha = alpha)) +
    geom_point() +
    scale_alpha_identity() +
    # ggtitle("IWS") +
    scale +
    xlab("Dim 3 (4.8%)") +
    ylab("Dim 4 (4.1%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.title = element_blank())

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_d34.png", width = 1500, height = 1200, res = 300)
iws_mfa_plot_d34
dev.off()


iws_mfa_plot_ori <- mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = oriCluster)) +
    geom_point() +
    # ggtitle("IWS") +
    xlab("Dim 1 (7.2%)") +
    ylab("Dim 2 (5.5%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_oriCluster.png", width = 1500, height = 1200, res = 300)
iws_mfa_plot_ori
dev.off()


iws_mfa_plot_ori_d34 <- mfa_df %>%
    ggplot(aes(x = PC3, y = PC4, color = oriCluster)) +
    geom_point() +
    # ggtitle("IWS") +
    xlab("Dim 3 (4.8%)") +
    ylab("Dim 4 (4.1%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_oriCluster_d34.png", width = 1500, height = 1200, res = 300)
iws_mfa_plot_ori_d34
dev.off()


IWS_sub_title <- ggdraw() + 
  draw_label("IWS - Sub-groups", 
             fontface = 'bold', x = 0.5, hjust = 0.5)

IWS_ori_title <- ggdraw() + 
  draw_label("IWS - Original Clusters", 
             fontface = 'bold', x = 0.5, hjust = 0.5)

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_panel.png", width = 1800, height = 1800, res = 300)
plot_grid(IWS_sub_title, 
    plot_grid(iws_mfa_plot, iws_mfa_plot_d34,  
        ncol = 2, rel_widths = c(1, 1.6)),
    IWS_ori_title,
    plot_grid(iws_mfa_plot_ori, iws_mfa_plot_ori_d34, ncol = 2, rel_widths = c(1, 1.6)),
    ncol = 1, rel_heights = c(0.1, 1))
dev.off()


## GESMD
gesmd_mfa_df <- tibble(Sample = gesmd_dataset_filt$ID,
    PC1 = GESMD_mfa$ind$coord[, 1],
    PC2 = GESMD_mfa$ind$coord[, 2],
    PC3 = GESMD_mfa$ind$coord[, 3],
    PC4 = GESMD_mfa$ind$coord[, 4]) %>%
    left_join(gesmd_dataset_filt, by = c("Sample" = "ID")) %>%
    mutate(sub_group = classifySamples(.),
    oriCluster = factor(clust_gesmd))

gesmd_mfa_plot <- gesmd_mfa_df %>%
    mutate(alpha = ifelse(sub_group %in% c("MDS-LB", "MDS-IB1", "MDS-IB2"), 0.3, 1)) %>%
    ggplot(aes(x = PC1, y = PC2, color = sub_group, alpha = alpha)) +
    geom_point() +
    scale_alpha_identity() +
    # ggtitle("GESMD") +
    scale +
    xlab("Dim 1 (6.7%)") +
    ylab("Dim 2 (5.1%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA.png", width = 1500, height = 1200, res = 300)
gesmd_mfa_plot
dev.off()

gesmd_mfa_plot_d34 <- gesmd_mfa_df %>%
    mutate(alpha = ifelse(sub_group %in% c("MDS-LB", "MDS-IB1", "MDS-IB2"), 0.3, 1)) %>%
    ggplot(aes(x = PC3, y = PC4, color = sub_group, alpha = alpha)) +
    geom_point() +
    scale_alpha_identity() +
    # ggtitle("GESMD") +
    scale +
    xlab("Dim 3 (4.9%)") +
    ylab("Dim 4 (4.3%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.title = element_blank())

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA_d34.png", width = 1500, height = 1200, res = 300)
gesmd_mfa_plot_d34
dev.off()


gesmd_mfa_ori_plot <- gesmd_mfa_df %>%
    ggplot(aes(x = PC1, y = PC2, color = oriCluster)) +
    geom_point() +
    # ggtitle("GESMD") +
    xlab("Dim 1 (6.7%)") +
    ylab("Dim 2 (5.1%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.position = "none")


png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA_oriCluster.png", width = 1500, height = 1200, res = 300)
gesmd_mfa_ori_plot
dev.off()

gesmd_mfa_ori_plot_d34 <- gesmd_mfa_df %>%
    ggplot(aes(x = PC3, y = PC4, color = oriCluster)) +
    geom_point() +
    # ggtitle("GESMD") +
    xlab("Dim 3 (4.9%)") +
    ylab("Dim 4 (4.3%)") +
    scale_color_discrete( name = "Original Cluster") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
    legend.title = element_blank())

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA_oriCluster_d34.png", width = 1500, height = 1200, res = 300)
gesmd_mfa_ori_plot_d34
dev.off()   

IWS_sub_title <- ggdraw() + 
  draw_label("IWS", 
             fontface = 'bold', x = 0.42, hjust = 0.5)


IWS_ori_title <- ggdraw() + 
  draw_label("IWS - Original Clusters", 
             fontface = 'bold', x = 0.5, hjust = 0.5)


GESMD_sub_title <- ggdraw() + 
  draw_label("GESMD", 
             fontface = 'bold', x = 0.42, hjust = 0.5)

GESMD_ori_title <- ggdraw() + 
  draw_label("GESMD - Original Clusters", 
             fontface = 'bold', x = 0.5, hjust = 0.5)


png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_GESMD_MFA_panel.png", width = 1800, height = 1800, res = 300)
plot_grid(IWS_sub_title, 
    plot_grid(iws_mfa_plot, iws_mfa_plot_d34,  
        ncol = 2, rel_widths = c(1, 1.6)),
    GESMD_sub_title,
    plot_grid(gesmd_mfa_plot, gesmd_mfa_plot_d34, ncol = 2, rel_widths = c(1, 1.6)),
    ncol = 1, rel_heights = c(0.1, 1))
dev.off()


png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_MFA_ori_panel.png", width = 1800, height = 900, res = 300)
plot_grid(IWS_ori_title, 
    plot_grid(iws_mfa_plot_ori, iws_mfa_plot_ori_d34,  
        ncol = 2, rel_widths = c(1, 1.6)),
    ncol = 1, rel_heights = c(0.1, 1))
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_MFA_ori_panel.png", width = 1800, height = 900, res = 300)
plot_grid(GESMD_ori_title,
    plot_grid(gesmd_mfa_ori_plot, gesmd_mfa_ori_plot_d34, ncol = 2, rel_widths = c(1, 1.6)),
    ncol = 1, rel_heights = c(0.1, 1))
dev.off()


## Features proportions cluster
features_clustIWS <- c("TET2bi", "EZH2",  "PHF6",  "STAG2", "Other Functions", "complex")
names(features_clustIWS) <- features_clustIWS
feat_tabs_clustIWS <- lapply(features_clustIWS, function(x) table(IWS_cluster[[x]], clust_comb1))
feat_tabs_clustIWS[["-7 or SETPBP1"]] <- table(IWS_cluster[["del7"]] == 1 | IWS_cluster[["SETBP1"]] == 1, clust_comb1)
feat_tabs_clustIWS[["del5q or SF3B1"]] <- table(IWS_cluster[["del5q"]] == 1 | IWS_cluster[["SF3B1"]] == 1, clust_comb1)
feat_tabs_props_clustIWS <- sapply(feat_tabs_clustIWS, function(x) prop.table(x, margin = 2)[2, ])*100

ph_IWS_oricluster_feat_prop <- pheatmap(mat = feat_tabs_props_clustIWS, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "IWS", legend = FALSE, silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_oricluster_feat_proportion.png", width = 1000, height = 1000, res = 300)
ph_IWS_oricluster_feat_prop
dev.off()

features_clustGESMD <- c("complex", "U2AF1", "STAG2")
names(features_clustGESMD) <- features_clustGESMD
feat_tabs_clustGESMD <- lapply(features_clustGESMD, function(x) table(GESMD_cluster[[x]], clust_gesmd))
feat_tabs_clustGESMD[["-7 or del5q"]] <- table(gesmd_mfa_df[["del7"]] == 1 | gesmd_mfa_df[["del5q"]] == 1, clust_gesmd)
feat_tabs_clustGESMD[["Hem. Diff. or SF3B1"]] <- table(gesmd_mfa_df[["Hematopoietic Differentiation"]] == 1 | gesmd_mfa_df[["SF3B1"]] == 1, clust_gesmd)
feat_tabs_clustGESMD[["TET2-bi or ZRSR2"]] <- table(gesmd_mfa_df[["TET2bi"]] == 1 | gesmd_mfa_df[["ZRSR2"]] == 1, clust_gesmd)
feat_tabs_clustGESMD[["EZH2 or SETBP1"]] <- table(gesmd_mfa_df[["EZH2"]] == 1 | gesmd_mfa_df[["SETBP1"]] == 1, clust_gesmd)
feat_tabs_props_clustGESMD <- sapply(feat_tabs_clustGESMD, function(x) prop.table(x, margin = 2)[2, ])*100

ph_GESMD_oricluster_feat_prop <- pheatmap(mat = feat_tabs_props_clustGESMD, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "GESMD", legend = FALSE, silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_oricluster_feat_proportion.png", width = 1200, height = 1300, res = 300)
ph_GESMD_oricluster_feat_prop
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_GESMD_oricluster_feat_proportion.png", width = 2200, height = 1000, res = 300)
plot_grid(ph_IWS_oricluster_feat_prop$gtable, 
    ph_GESMD_oricluster_feat_prop$gtable, 
    ncol = 2)
dev.off()

## Features proportions
features <- c("TET2bi", "EZH2",  "del7", "STAG2", "del5q", "SF3B1", "complex")
names(features) <- features
feat_tabs <- lapply(features, function(x) table(mfa_df[[x]], mfa_df$sub_group))
feat_tabs_props <- sapply(feat_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

colnames(feat_tabs_props) <- c("TET2-bi", "EZH2", "-7", "STAG2", "del5q", "SF3B1", "Complex Kar.")

ph_IWS_feat_prop <- pheatmap(mat = feat_tabs_props, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "IWS", legend = FALSE, silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_feat_proportion.png", width = 1000, height = 1000, res = 300)
ph_IWS_feat_prop
dev.off()

gesmd_feat_tabs <- lapply(features, function(x) table(gesmd_mfa_df[[x]], gesmd_mfa_df$sub_group))
gesmd_feat_tabs_props <- sapply(gesmd_feat_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
colnames(gesmd_feat_tabs_props) <- c("TET2-bi", "EZH2", "-7", "STAG2", "del5q", "SF3B1", "Complex Kar.")

ph_GESMD_feat_prop <- pheatmap(mat = gesmd_feat_tabs_props, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, main = "GESMD", legend = FALSE, silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_feat_proportion.png", width = 1000, height = 1000, res = 300)
ph_GESMD_feat_prop
dev.off()


png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_GESMD_feat_proportion.png", width = 2200, height = 1000, res = 300)
plot_grid(ph_IWS_feat_prop$gtable, ph_GESMD_feat_prop$gtable, ncol = 2)
dev.off()


## Selected components
png("figures/GESMD_IWS_clustering/subgroup_definition/joint_components_inertia.png", width = 400, height = 300)
tibble(Dimension = seq_len(40), IWS = IWS_mfa$global.pca$eig[1:40, 2], GESMD = GESMD_mfa$global.pca$eig[1:40, 2]) %>%
    pivot_longer(cols = c("IWS", "GESMD"), names_to = "Dataset", values_to = "Inertia") %>%
    mutate(Selection = ifelse((Dimension <= 9 & Dataset == "IWS") | (Dimension <= 6 & Dataset == "GESMD"), "chosen", "rest"),
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
    mutate(color = ifelse((N_clusters == 10 & Dataset == "IWS") | (N_clusters == 10 & Dataset == "GESMD"), "selected", "rest")) %>%
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

ph_IWS <- pheatmap(mat = IWS_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0,
    main = "IWS", silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_cluster_mapping.png", width = 1800, height = 1500, res = 300)
ph_IWS
dev.off()

GESMD_tab <- table(gesmd_mfa_df$sub_group, clust_gesmd)
GESMD_tab <- GESMD_tab[levels(gesmd_mfa_df$sub_group), ]


ph_gesmd <- pheatmap(mat = GESMD_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0,
    main = "GESMD", silent = TRUE)

png("figures/GESMD_IWS_clustering/subgroup_definition/GESMD_cluster_mapping.png", width = 1800, height = 1500, res = 300)
ph_gesmd
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_definition/IWS_GESMD_cluster_mapping.png", width = 1800, height = 2000, res = 300)
plot_grid(ph_IWS$gtable, ph_gesmd$gtable, ncol = 1)
dev.off()