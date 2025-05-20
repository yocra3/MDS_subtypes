#' ---------------------------
#'
#' Purpose of script:
#'
#' Explore clusters obtained with mutations with umap
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
library(survival)
library(survminer)
library(pheatmap)
library(cowplot)
library(umap)
library(MASS)

load( "results/clustering/umap_mutations.Rdata")
load("results/clustering/umap_mutations_clustering.Rdata")
clinical_clust <- clinical_clust %>% 
    mutate(clust = factor(cluster),
    complex = factor(complex, levels = c("non-complex", "complex")),
    mutation0 = as.numeric(N_mutations == 0),
    mutation1 = as.numeric(N_mutations == 1),
    mutation2_3 = as.numeric(N_mutations %in% c(2, 3)),
    mutation_multi = as.numeric(N_mutations > 3))  %>%
    filter(!is.na(clust))
cluster_levs <- levels(clinical_clust$clust)
mds_types <- unique(clinical_clust$consensus)
cluster_map <- select(clinical_clust, consensus, clust) %>%
    distinct()

#' ---------------------------
# Explore clustering features 
#' ---------------------------

## Sex
#' ---------------------------
table(clinical_clust$clust, clinical_clust$SEX)
sex_prop <- prop.table(table(clinical_clust$clust, clinical_clust$SEX), margin = 1)
sex_prop_df <- tibble(cluster = rownames(sex_prop), propF = sex_prop[, 1])
chisq.test(table(clinical_clust$cluster, clinical_clust$SEX))

sex_test <- lapply(levels(clinical_clust$clust), function(x) {
    test <- fisher.test(table(clinical_clust$clust == x, clinical_clust$SEX))
    c(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
} ) %>%
    Reduce(rbind, .) %>%
    as_tibble() %>%
    mutate(cluster = as.character(cluster),
            pval = as.numeric(pval),
            OR = as.numeric(OR))

png("figures/umap_subtypes/clust_Sex.png")
left_join(sex_prop_df, sex_test, by = "cluster") %>%
    mutate(Significant = ifelse(pval < 0.05/length(cluster_levs), "Significant", "Non-significant"),
        clust = factor(cluster)) %>%
    left_join(cluster_map, by = "clust") %>%
    ggplot(aes(y = propF*100, x = cluster, fill = Significant)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        geom_hline(yintercept = mean(clinical_clust$SEX == "F")*100,
            linetype = "dashed") +
            ylab("Proportion of Females") +
            facet_wrap(~ consensus, scales = "free_x")
dev.off()


## Age
#' ---------------------------
tapply(clinical_clust$AGE, clinical_clust$clust, summary)

png("figures/umap_subtypes/clust_age.png")
ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_clust$AGE, na.rm = TRUE),
        linetype = "dashed") +
            facet_wrap(~ consensus, scales = "free_x")
dev.off()


## IPSSM
#' ---------------------------
table(clinical_clust$clust, clinical_clust$IPSSM)
prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1)
chisq.test(table(clinical_clust$clust, clinical_clust$IPSSM))

png("figures/umap_subtypes/clust_ipssm.png")
table(clinical_clust$clust, clinical_clust$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/umap_subtypes/clust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()


## Kariotipos
#' ---------------------------
kar_events <- c("delY", "del11q", "del5q", "del12p",
                "del20q", "del7q", "plus8", "plus19", "del7", "complex")
names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(clinical_clust[[x]], clinical_clust$clust))

png("figures/umap_subtypes/clust_karevents_prop.png")
kar_tabs_prop <- sapply(kar_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = kar_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/clust_karevents_prop2.png")
kar_tabs_prop2 <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = kar_tabs_prop2, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/lowblasts_karevents_prop2.png")
kar_tabs_lowblasts <- kar_tabs_prop2[paste0("Low blasts_", 1:4), ]
pheatmap(mat = kar_tabs_lowblasts, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/ib1_karevents_prop2.png")
kar_tabs_ib1 <- kar_tabs_prop2[paste0("MDS-IB1_", 1:3), ]
pheatmap(mat = kar_tabs_ib1, display_numbers = TRUE)
dev.off()



## Mutaciones
#' ---------------------------
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

mutations <- colnames(clinical_clust)[colnames(clinical_clust) %in% colnames(mutation)]
sel_muts <- mutations[mutations != "ID"]
names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(clinical_clust[[x]], clinical_clust$clust))


png("figures/umap_subtypes/clust_mutation_prop.png", width = 2000, height = 800)
mut_tabs_prop <- sapply(mut_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = mut_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/clust_mutation_prop2.png", width = 2000, height = 800)
mut_tabs_prop2 <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = mut_tabs_prop2, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/lowblasts_mutation_prop2.png")
mut_tabs_lowblasts <- mut_tabs_prop2[paste0("Low blasts_", 1:4), ]
pheatmap(mat = mut_tabs_lowblasts[, apply(mut_tabs_lowblasts, 2, max) > 10], display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/ib1_mutation_prop2.png")
mut_tabs_ib1 <- mut_tabs_prop2[paste0("MDS-IB1_", 1:3), ]
pheatmap(mat = mut_tabs_ib1[, apply(mut_tabs_ib1, 2, max) > 10], display_numbers = TRUE)
dev.off()


## Variables clinicas
#' ---------------------------
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clinical_clust, aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(~ consensus, scales = "free_x")
})
png("figures/umap_subtypes/clust_clinical_vars.png", width = 2000, height = 1000)
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()

clin_plots_lowblasts <- lapply(clin_vars, function(x){
    clinical_clust %>%
        filter(consensus == "Low blasts") %>%
        ggplot(aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw() 
})

png("figures/umap_subtypes/lowblasts_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots_lowblasts, nrow = 3)
dev.off()

clin_plots_ib1 <- lapply(clin_vars, function(x){
    clinical_clust %>%
        filter(consensus == "MDS-IB1") %>%
        ggplot(aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw() 
})

png("figures/umap_subtypes/ib1_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots_ib1, nrow = 3)
dev.off()


#' ---------------------------
# Effect of clusters on survival 
#' ---------------------------

## Comparison with survival
#' ---------------------------
cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_clust)

clinical_clust$clust <- relevel(clinical_clust$clust, "6") ## Lowest scores

png("figures/umap_subtypes/clust_survival.png", height = 800, width = 2000)
a <- lapply(mds_types, function(mds){
p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust, subset = consensus == mds) %>%
    ggsurvplot(data = clinical_clust, title = mds)
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

png("figures/umap_subtypes/lowblasts_survival.png", height = 600, width = 1000)
p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust, subset = consensus == "Low blasts") %>%
    ggsurvplot(data = clinical_clust, title = "Low blasts")
    p$plot + theme_bw(base_size = 20) 
dev.off()


png("figures/umap_subtypes/ib1_survival.png", height = 600, width = 1000)
p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust, subset = consensus == "MDS-IB1") %>%
    ggsurvplot(data = clinical_clust, title = "IB1")
    p$plot + theme_bw(base_size = 20) 
dev.off()



lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ cluster + AGE + SEX, clinical_clust,subset = consensus == mds )

})

lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ cluster + AGE + SEX + IPSSM_SCORE, clinical_clust,subset = consensus == mds )

})

## Comparison with transformation to leukemia
#' ---------------------------
lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ cluster + AGE + SEX, clinical_clust,subset = consensus == mds )

})

lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ cluster + AGE + SEX + IPSSM_SCORE, clinical_clust,subset = consensus == mds )

})

## Comparison with Leukemia Free Survival
#' ---------------------------
lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(LFS_YEARS,LFS_STATUS) ~ cluster + AGE + SEX, clinical_clust,subset = consensus == mds )

})

## Comparison with treatment
#' ---------------------------
lapply(mds_types, function(mds){

    cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ cluster*hma + AGE + SEX, clinical_clust,subset = consensus == mds )

})


## LOW blasts abstract 
#' ---------------------------
set.seed(27)
cluster_low <- lapply("Low blasts", function(mds){
    k <- 4
    umap_res <- umap_types[[mds]]
    # umap_dist <- dist(umap_res$layout)
    # hclust <- hclust(umap_dist)

    # clust_class <- cutree(hclust, k = k)
    km <- kmeans(umap_res$layout, centers = k, nstart = 1000)
    clust_class <- km$clust
    paste(mds, clust_class, sep = "_")
})
table(cluster_low)

umap_low <- as_tibble(umap_types[["Low blasts"]]$layout) %>%
        mutate(Clusters = factor(unlist(cluster_low)))
levels(umap_low$Clusters) <- c("Mutational", "AMML-like", "Basal", "Aneuploidies")
umap_low$Clusters <- factor(umap_low$Clusters, levels = c("Basal", "Aneuploidies", "Mutational", "AMML-like"))

low_umap_plot <- ggplot(umap_low, aes(x = V1, y = V2, col = Clusters)) +
        geom_point(alpha = 0.5) +
        theme_bw()
png("figures/umap_subtypes/lowblasts_umap.png", height = 600, width = 600, res = 150)
low_umap_plot + ggtitle("A")
dev.off()



clinical_low <- clinical_clust %>%
        filter(consensus == "Low blasts" & !is.na(cluster)) %>%
        mutate(cluster = factor(cluster),
               cluster = fct_recode(cluster, 
                        "Basal" = "Low blasts_1",
                        "Aneuploidies"= "Low blasts_2",
                        "Mutational" = "Low blasts_3",
                        "AMML-like" = "Low blasts_4")
        )
save(clinical_low, file = "results/clustering/clinical_lowblasts.Rdata")

clinical_low$logPLT <- log10(clinical_low$PLT)
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "logPLT")
names(clin_vars) <- clin_vars
low_clin_plots <- lapply(clin_vars, function(var){
    ggplot(clinical_low, aes(x = .data[[var]])) +
        geom_histogram() +
        facet_grid(cluster ~ .) +
        theme_bw()
})
png("figures/umap_subtypes/lowblasts_clinical.png", width = 1000, height = 1000)
plot_grid(plotlist = low_clin_plots, nrow = 3)
dev.off()

lapply(clin_vars[clin_vars != "logPLT"], function(var){
    summary(glm(formula (paste(var, " ~ cluster")), clinical_low, family = "poisson"))
})


summary(lm(logPLT ~ cluster, clinical_low)) 


clinvars_low <- clinical_low %>%
    select(BM_BLAST, PB_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
    as.matrix()

pcclin_low <- prcomp(clinvars_low, scale = TRUE)
png("figures/umap_subtypes/lowblasts_clinical_pc.png", width = 1000)
p1 <- data.frame(PC1 = pcclin_low$x[, 1], PC2 = pcclin_low$x[, 2], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = PC1, y = PC2, col = cluster)) +
        geom_point() +
        theme_bw()
p2 <- data.frame(PC3 = pcclin_low$x[, 3], PC4 = pcclin_low$x[, 4], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = PC3, y = PC4, col = cluster)) +
        geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()

tapply(pcclin_low$x[, 1], clinical_low$cluster, summary)
tapply(pcclin_low$x[, 2], clinical_low$cluster, summary)

summary(lm(pcclin_low$x[, 1] ~ clinical_low$cluster))
summary(lm(pcclin_low$x[, 2] ~ clinical_low$cluster))

clin_dist <- dist(clinvars_low, method = "euclidean")
clin_dist_pca <- isoMDS(clin_dist, k = 5)
cor(clin_dist_pca$points, pcclin_low$x)

png("figures/umap_subtypes/lowblasts_clinical_distpc.png", width = 1000)
p1 <- data.frame(MDS1 = clin_dist_pca$points[, 1], MDS2 = clin_dist_pca$points[, 2], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = MDS1, y = MDS2, col = cluster)) +
        geom_point() +
        theme_bw()
p2 <- data.frame(MDS3 = clin_dist_pca$points[, 3], MDS4 = clin_dist_pca$points[, 4], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = MDS3, y = MDS4, col = cluster)) +
        geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()


kar_mat <- clinical_low[, kar_events] %>% as.matrix()
kar_mat  <- kar_mat [, colMeans(kar_mat ) > 0.01] ## Remove low frequency karyotypes
kar_mat <- (kar_mat - 1) * -1
kar_pca <- ca(kar_mat)

png("figures/umap_subtypes/lowblasts_kar_distpc.png", width = 1000)
p1 <- data.frame(CA1 = kar_pca$rowcoord[, 1], CA2 = kar_pca$rowcoord[, 2], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = CA1, y = CA2, col = cluster)) +
   #     geom_point() +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        theme_bw()
p2 <- data.frame(CA3 = kar_pca$rowcoord[, 3], CA4 = kar_pca$rowcoord[, 4], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = CA3, y = CA4, col = cluster)) +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        #geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()

mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
mutations <- colnames(clinical_low)[colnames(clinical_low) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

mut_mat <- clinical_low[, mut_vars] %>% as.matrix()
mut_mat <- mut_mat[, colMeans(mut_mat) > 0.05]
mut_mat <- (mut_mat - 1) * -1
mut_pca <- ca(mut_mat)

png("figures/umap_subtypes/lowblasts_mut_distpc.png", width = 1000)
p1 <- data.frame(CA1 = mut_pca$rowcoord[, 1], CA2 = mut_pca$rowcoord[, 2], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = CA1, y = CA2, col = cluster)) +
     #   geom_point() +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        theme_bw()
p2 <- data.frame(CA3 = mut_pca$rowcoord[, 3], CA4 = mut_pca$rowcoord[, 4], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = CA3, y = CA4, col = cluster)) +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        #geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()

dist_mut_pc <- dist(mut_pca$rowcoord[, 1:4]*matrix(mut_pca$sv[1:4], ncol = 4, nrow = nrow(mut_pca$rowcoord), byrow = TRUE),
     method = "euclidean")
dist_mut_hc <- hclust(dist_mut_pc)
clust_mut <- cutree(dist_mut_hc, k = 10)
table(clust_mut, clinical_low$cluster)


png("figures/umap_subtypes/lowblasts_mut_distpc2.png", width = 1000)
p1 <- data.frame(CA1 = mut_pca$rowcoord[, 1], CA2 = mut_pca$rowcoord[, 2], cluster = factor(clust_mut)) %>%
    ggplot(aes(x = CA1, y = CA2, col = cluster)) +
     #   geom_point() +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        theme_bw()
p2 <- data.frame(CA3 = mut_pca$rowcoord[, 3], CA4 = mut_pca$rowcoord[, 4], cluster = factor(clust_mut)) %>%
    ggplot(aes(x = CA3, y = CA4, col = cluster)) +
        geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
        #geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()


comb_pc_mat <- cbind(pcclin_low$x[, 1:4], kar_pca$rowcoord[, 1:3], mut_pca$rowcoord[, 1:16])
comb_pc_mat <- cbind(pcclin_low$x[, 1:4], kar_pca$rowcoord[, 1:3], mut_pca$rowcoord[, 1:5])

comb_pc <- prcomp(comb_pc_mat)

png("figures/umap_subtypes/lowblasts_combined_pc.png", width = 1000)
p1 <- data.frame(PC1 = comb_pc$x[, 1], PC2 = comb_pc$x[, 2], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = PC1, y = PC2, col = cluster)) +
        geom_point() +
        theme_bw()
p2 <- data.frame(PC3 = comb_pc$x[, 3], PC4 = comb_pc$x[, 4], cluster = clinical_low$cluster) %>%
    ggplot(aes(x = PC3, y = PC4, col = cluster)) +
        geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()

lapply(2:12, function(i){
    dist_comb_pc <- dist(comb_pc$x[, 1:i], method = "euclidean")
dist_comb_hc <- hclust(dist_comb_pc)
clust_comb <- cutree(dist_comb_hc, k = 8)
table(clust_comb, clinical_low$cluster)
})
clust_comb <- kmeans(comb_pc$x[, 1:8], centers = 4, nstart = 1000)$cluster
dist_comb_pc <- dist(comb_pc$x[, 1:8], method = "euclidean")
table(clust_comb, clinical_low$cluster)

png("figures/umap_subtypes/lowblasts_combined_pc2.png", width = 1000)
p1 <- data.frame(PC1 = comb_pc$x[, 1], PC2 = comb_pc$x[, 2], cluster = factor(clust_comb6)) %>%
    ggplot(aes(x = PC1, y = PC2, col = cluster)) +
        geom_point() +
        theme_bw()
p2 <- data.frame(PC3 = comb_pc$x[, 3], PC4 = comb_pc$x[, 4], cluster = factor(clust_comb6)) %>%
    ggplot(aes(x = PC3, y = PC4, col = cluster)) +
        geom_point() +
        theme_bw()
plot_grid(p1, p2, ncol = 2)
dev.off()

library(cluster)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

sil_scores <- sapply(2:40, silhouette_score,
      mat = comb_pc$x[, 1:8], data_dist = dist_comb_pc)


png("figures/lowblast_hc_cluster/comb_kmeans_silhouette.png")
plot(2:40, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()
clust_comb6 <- kmeans(comb_pc$x[, 1:8], centers = 6, nstart = 1000)$cluster
clust_comb14 <- kmeans(comb_pc$x[, 1:8], centers = 14, nstart = 1000)$cluster
table(clust_comb6, clinical_low$cluster)
table(clust_comb14, clinical_low$cluster)

library(rpart)

clinical_low_tree <- dplyr::select(clinical_low, 
    all_of(c(clin_vars, "AGE", "SEX", kar_events, mut_vars))) %>%
    mutate(clust = clust_comb6)

tree_mod <- rpart(clust ~ . , data = clinical_low_tree, method = "class")

predictions <- predict(tree_mod, clinical_low_tree, type = "class")
table(predictions, clust_comb6)

clinical_low_tree$clust <- predictions
tree_mod <- rpart(clust ~ . , data = clinical_low_tree, method = "class")
predictions2 <- predict(tree_mod, clinical_low_tree, type = "class")
table(predictions, predictions2)

clinical_low_tree$clust <- predictions2
save(clinical_low_tree, file = "results/clustering/clinical_lowblasts_tree.Rdata")

clinical_low$clust <- factor(predictions2)
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "logPLT")
names(clin_vars) <- clin_vars
low_clin_plots2 <- lapply(clin_vars, function(var){
    ggplot(clinical_low, aes(x = .data[[var]])) +
        geom_histogram() +
        facet_grid(clust ~ .) +
        theme_bw()
})
png("figures/umap_subtypes/lowblasts_clinical2.png", width = 1000, height = 1000)
plot_grid(plotlist = low_clin_plots2, nrow = 3)
dev.off()

lapply(clin_vars[clin_vars != "logPLT"], function(var){
    summary(glm(formula (paste(var, " ~ clust")), clinical_low, family = "poisson"))
})
summary(lm(logPLT ~ clust, clinical_low)) 

png("figures/umap_subtypes/lowblast_clust_ipssm.png")
table(clinical_low$clust, clinical_low$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/umap_subtypes/lowblast_clust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_low$clust, clinical_low$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()


cox_surv_raw <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_low )
cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX + IPSSM_SCORE, clinical_low )

surv_low <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/umap_subtypes/lowblasts_survival_new.png")
surv_low$plot 
dev.off()

cox_aml_raw <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust + AGE + SEX , clinical_low )
cox_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust + AGE + SEX + IPSSM_SCORE, clinical_low )
aml_low <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ clust, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/umap_subtypes/lowblasts_amlt_new.png")
aml_low$plot
dev.off()


clinical_low$complex <- factor(clinical_low$complex, levels = c("non-complex", "complex"))
names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(clinical_low[[x]], clinical_low$clust))
kar_tabs_lowblasts <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/umap_subtypes/lowblasts_karevents_prop2_new.png")
pheatmap(mat = kar_tabs_lowblasts, display_numbers = TRUE)
dev.off()

sel_muts <- colnames(mut_mat)
names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(clinical_low[[x]], clinical_low$clust))

mut_tabs_lowblasts <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
png("figures/umap_subtypes/lowblasts_mutevents_prop2_new.png")
pheatmap(mat = mut_tabs_lowblasts, display_numbers = TRUE)
dev.off()

png("figures/umap_subtypes/lowblasts_classification_tree.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()


dist_comb_pc <- dist(comb_pc$x[, 1:9], method = "euclidean")
dist_comb_hc <- hclust(dist_comb_pc)
clust_comb <- cutree(dist_comb_hc, k = 8)
table(clust_comb, clinical_low$cluster)

cor(pcclin_low$x, clinvars_low)


clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    mutate(clust = relevel(clust, 4)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    mutate(clust = relevel(clust, 2)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ STAG2*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
     mutate(clust = relevel(clust, 6)) %>%
     coxph(Surv(OS_YEARS,OS_STATUS) ~ EZH2*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(2, 3, 5)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ZRSR2*clust + AGE + SEX + IPSSR_SCORE, data = .) 


clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*clust + AGE + SEX + IPSSR_SCORE, data = .) 

############# Old code
















surv_low$plot + theme_bw(base_size = 15) + ggtitle("C") +
    geom_label(aes(x = 12, y = 1, label = "P = 1.5e-4"), size = 5, color = "black")
dev.off()

clin_plots_lowblasts <- clinical_low %>%
        ggplot(aes(y = MONOCYTES, x = cluster)) +
    geom_boxplot() +
    theme_bw() 

cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ cluster + AGE + SEX, clinical_low )

surv_low <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ cluster, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/umap_subtypes/lowblasts_survival2.png", height = 600, width = 800, res = 150)
surv_low$plot + theme_bw(base_size = 15) + ggtitle("C") +
    geom_label(aes(x = 12, y = 1, label = "P = 1.5e-4"), size = 5, color = "black")
dev.off()

cox_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ cluster + AGE + SEX, clinical_low )

aml_low <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ cluster, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/umap_subtypes/lowblasts_amlt.png",  height = 600, width = 800, res = 150)
aml_low$plot + theme_bw(base_size = 15) +
 ylab("AML transformation") + ggtitle("D") +
 geom_label(aes(x = 12, y = 1, label = "P = 1.5e-4"), size = 5, color = "black")
dev.off()



## Kariotipos
#' ---------------------------
kar_events <- c("plus8", "del7", "complex")
clinical_low$complex <- factor(clinical_low$complex, levels = c("non-complex", "complex"))
names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(clinical_low[[x]], clinical_low$cluster))
kar_tabs_lowblasts <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/umap_subtypes/lowblasts_karevents_prop2.png")
pheatmap(mat = kar_tabs_lowblasts, display_numbers = TRUE)
dev.off()

## Mutaciones
#' ---------------------------
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

mutations <- colnames(clinical_clust)[colnames(clinical_clust) %in% colnames(mutation)]
sel_muts <- c("TET2", "ASXL1", "DNMT3A", "RUNX1", "SRSF2")
names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(clinical_low[[x]], clinical_low$cluster))

mut_tabs_lowblasts <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
comb_tabs <- cbind(kar_tabs_lowblasts, mut_tabs_lowblasts)

kar_mut_plot <- pheatmap(mat = round(comb_tabs, 0), display_numbers = TRUE, legend_labels = "Legend Title")
heatmap_grob <- kar_mut_plot$gtable

library(grid)
library(gridExtra)
title_grob <- textGrob("B", gp = gpar(fontsize = 14, fontface = "bold"), 
    x = 0.1, hjust = 0, y = 0.9)

legend_grob <- textGrob("Cluster Prop.", rot = -90, gp = gpar(fontsize = 14), 
    x = 0.3)


png("figures/umap_subtypes/lowblasts_karmut_prop.png", height = 600, width = 800, res = 150)
grid.newpage()
grid.arrange(title_grob, heatmap_grob, legend_grob,
              ncol = 3, widths = c(1, 20, 1)
     )
dev.off()

png("figures/umap_subtypes/abstract_panel.png", width = 600, height = 800, res = 300)
plot_grid(low_umap_plot, heatmap_grob, 
    surv_low$plot + theme_bw(base_size = 3),
    aml_low$plot + theme_bw(base_size = 3) + ylab("AML transformation"), 
    nrow = 4, labels = "AUTO")
dev.off()




low_subtypes <- c("Basal", "Aneuploidies", "AMML", "Mutational")
sel_genes <- c("RUNX1", "SRSF2", "DNMT3A", "ASXL1", "U2AF1", "TET2", "ZRSR2", "STAG2", "EZH2")
mut_hr_clust <- lapply(sel_genes, function(gene){
    lapply(low_subtypes, function(x){
        mod <- coxph(formula(paste("Surv(OS_YEARS,OS_STATUS) ~", gene, " + AGE + IPSSR_SCORE")), 
            clinical_low, subset = cluster == x)
        coef <- summary(mod)$coefficients
        hr <- exp(coef[1, 1])
        pval <- coef[1, 5]
        freq <- mean(clinical_low[[gene]][clinical_low$cluster == x], na.rm = TRUE)
        data.frame(HR = hr, Pvalue = pval, Cluster = x, Gene = gene, Freq = freq )
    }) %>% Reduce(f = rbind) %>%
        as_tibble()
}) %>% Reduce(f = rbind) %>%
    mutate(Cluster = factor(Cluster, levels = low_subtypes))



## Modify one value due to low frequency
mut_hr_clust <- mut_hr_clust %>%
    mutate(HR = ifelse(HR < 0.59, 0.59, HR))

## add colors
mut_hr_col <- mutate(mut_hr_clust, 
    Sig = ifelse(Pvalue < 0.05, "Signif", "Not-signif"),
    Direction = ifelse(HR > 1, "Risk", "Protective"), 
    Color = ifelse(Freq < 0.1, "Undetermined", paste(Direction, Sig)),
    Color = factor(Color, levels = c("Protective Signif", "Protective Not-signif",
        "Risk Signif", "Risk Not-signif", "Undetermined")) )


mut_hr_plot <- mut_hr_col %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        facet_wrap(~ Gene) +
        scale_fill_manual(name = "", values = c("green", "darkgreen", "red", "darkred", "grey"))


png("figures/umap_subtypes/mutations_risk.png", width = 3000, height = 1000, res = 300)
mut_hr_plot
dev.off()

clinical_low %>%
    mutate(Cluster = ifelse(cluster == "Aneuploidies", "Aneuploidies", "Other")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*Cluster + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
filter(cluster != "Basal") %>%
   mutate(Cluster = ifelse(cluster != "Mutational", "Other", cluster)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ TET2*Cluster + AGE + SEX  + IPSSR_SCORE, data = .) 

clinical_low %>%
    filter(cluster %in% c("Aneuploidies", "Mutational")) %>%
   mutate(Cluster = relevel(cluster, "Mutational"),
            Cluster = droplevels(Cluster)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ STAG2*Cluster + SEX + + AGE + IPSSR_SCORE, data = .) 

clinical_low %>%
    filter(!cluster %in% c("AMML")) %>%
   mutate(Cluster = ifelse(cluster != "Aneuploidies", "Other", cluster)) %>% 
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*Cluster + SEX + AGE + IPSSR_SCORE, data = .) 

coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*cluster + AGE + IPSSM_SCORE, data = clinical_low) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*cluster + SEX + AGE + IPSSR_SCORE, data = clinical_low) 
