#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define clusters of patients using FactoMineR
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.1 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(FactoMineR)
library(survival)
library(survminer)
library(missMDA)

load("results/preprocess/clinical_preproc.Rdata")

## Use only clinical variables for clustering
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clust_dataset <- clinical_all %>%
    select(ID, clin_vars) %>%
        mutate(logPLT = log(PLT),
            rat_BM_HB = BM_BLAST/HB,
            rat_WBC_ANC = WBC/(ANC + 0.1),
            rat_WBC_HB = WBC/HB,
            rat_WBC_PLT = WBC/logPLT) 


## Impute values
clin_impute_l <- imputePCA(clust_dataset[, -1])
clin_impute <- clin_impute_l$completeObs
colnames(clin_impute) <- colnames(clust_dataset)[-1]
rownames(clin_impute) <- clust_dataset$ID
## Run PCA with FactoMiner
pca_clin <- PCA(clin_impute )

## Make clusters with Hierarchical Clustering on Principal Components (HCPC) (original results)
pca_clust_ori <- HCPC(pca_clin,  nb.clust = 14)

png("figures/clin_clust_inertia.png")
plot(pca_clust$call$t$inert.gain[1:30])
dev.off()

clin_clust_ori <- pca_clust_ori$data.clust$clust

## Make clusters with HCPC - without k-means
pca_clust14 <- HCPC(pca_clin,  nb.clust = 14, consol = FALSE)
clin_clust14 <- pca_clust14$data.clust$clust
table(clin_clust_ori, clin_clust14)

lapply(2:13, function(nb){
    pca_clust <- HCPC(pca_clin,  nb.clust = nb, consol = FALSE)
    clin_clust <- pca_clust$data.clust$clust
    table(clin_clust, clin_clust14)
})
save(clin_clust_ori, pca_clust_ori, pca_clin, clin_impute_l, file = "results/clustering/clinical_clustering_ori.Rdata")
save(pca_clust14, clin_clust14,  file = "results/clustering/clinical_clustering_14.Rdata")



pca_clust6 <- HCPC(pca_clin,  nb.clust = 6, consol = FALSE)
clin_clust6 <- pca_clust6$data.clust$clust
save(pca_clust6, clin_clust6,  file = "results/clustering/clinical_clustering_6.Rdata")


pca <- pca_clin$ind$coord %>%
    as_tibble() %>%
    mutate(clust = factor(clin_clust6, levels = 1:6))

pdf("figures/pca_clust6.pdf")
ggplot(pca, aes(x = `Dim.1`, y = `Dim.2`, color = clust)) +
        geom_point(alpha = 0.4) +
        theme_bw()
ggplot(pca, aes(x = `Dim.3`, y = `Dim.4`, color = clust)) +
        geom_point(alpha = 0.4) +
        theme_bw()
ggplot(pca, aes(x = `Dim.4`, y = `Dim.5`, color = clust)) +
        geom_point(alpha = 0.4) +
        theme_bw()
dev.off()

