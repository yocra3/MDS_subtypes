#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define clusters of patients using UMAP
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.4 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(missMDA)
library(umap)
library(cluster)

load("results/preprocess/clinical_preproc.Rdata")
## Select variables for clustering
#' - Demographic: age 
#' - Clinical: Blood cell proportions
#' - Karyotype events: number of aberrations 
demo_vars <- c("AGE", "Sex")
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("monoy", "del11q", "del5q", "del12p",
                "del20q", "del7q", "triso8", "triso19", "mono7", 
                "i17q", "chr3ab")


clust_dataset <- clinical_all %>%
    mutate(Sex = ifelse(SEX == "F", 1, 0),
        PLT = log(PLT)) %>% ## Convert PLT to log
    select(c(Sex, AGE, clin_vars, kar_events, "N_aberrations"))
impute_l <- imputePCA(clust_dataset)
impute_mat <- impute_l$completeObs

cont_mat <- impute_mat[, c("AGE", clin_vars, "N_aberrations")]
cont_mat_scale <- scale(cont_mat)

clust_mat <- cbind(cont_mat_scale, impute_mat[, c("Sex", kar_events)])
rownames(clust_mat) <- clinical_all$ID

silhouette_score_tree <- function(data_dist, clusth, k) {
  km <- cutree(clusth, k = k)
  ss <- silhouette(km, data_dist)
  mean(ss[, 3])
}

getClusters <- function(neighbors, comps = 2, dist){
     
     set.seed(27)
     umap_res <- umap(clust_mat, n_neighbors = neighbors, n_components = comps, min_dist = dist)

     png(paste0("figures/UMAPs/umap_", comps, "comp_", neighbors, "neigh_", dist, "d.png"))
     plot(umap_res$layout)
     dev.off()

     umap_dist <- dist(umap_res$layout)
     umap_hclust <- hclust(umap_dist)

     sil_scores_tree <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist, clusth = umap_hclust)
     
     png(paste0("figures/UMAPs/umap_", comps, "comp_", neighbors, "neigh_", dist, "d_silhouettetree.png"))
     plot(2:40, sil_scores_tree, type="b", pch = 19, frame = FALSE,
          xlab = "Number of clusters K",
          ylab = "Silhouette Score",
          main = "Silhouette Analysis")
     dev.off()
     umap_hclust
}

## Run UMAP in bucle
neighbors <- c(15, 30, 50, 100)
components <- c(2, 6, 8)
dists <- c(0.001, 0.01, 0.1)

combinations <- expand.grid(neighbors, components, dists)

clust_list <- lapply(seq_len(nrow(combinations)), function(i){
     print(combinations[i, ])
     getClusters(combinations[i, 1], combinations[i, 2], combinations[i, 3])
})
save(combinations, clust_list, file = "results/clustering/umap_cluster_combs.Rdata")

## Look at plots and select most interesting
sel_combs <- list(c(1, 23), c(2, 15), c(3, 15), c(5, 17), c(6, 19),
     c(7, 17), c(9, 20), c(10, 19), c(11, 14), c(12, 13), c(13, 17),
     c(14, 18), c(15, 18), c(17, 21), c(18, 20), c(19, 15), c(21, 22),
     c(22, 18), c(23, 20), c(25, 22), c(26, 18), c(29, 15), c(30, 20), 
     c(33, 16), c(34, 18), c(35, 18))

umap_class_mat <- sapply(sel_combs, function(vec){
     cutree(clust_list[[vec[1]]], k = vec[2])
})

umap_dist_mat <- apply(umap_class_mat, 1, function(x) {
    apply(umap_class_mat, 1, function(y) sum(x != y))
})
meta_umap <- hclust(as.dist(umap_dist_mat))
meta_class <- cutree(meta_umap, h = 6) ## Samples in the same cluster 80% of time

load( "results/tsne_perp10.Rdata")

tb <- tibble(meta = meta_class, tsne10 = clust10_tree)
n_meta <- group_by(tb, meta) %>%
     summarize(n = n()) %>%
     filter(n > 20)

tb_filt <- mutate(tb, meta2 = ifelse(!meta %in% n_meta$meta, NA, meta))
table(tb_filt$meta2, tb_filt$tsne10, useNA = "ifany")
round(
     prop.table(
          table(tb_filt$meta2, tb_filt$tsne10, useNA = "ifany"),
      margin = 2) * 100,
1 )
   
meta50 <- filter(n_meta, n >= 50)$meta


clinical_tree <- select(clinical_all, c(clin_vars, "AGE", "SEX", kar_events, "N_aberrations")) %>%
    mutate(meta = meta_class,
          tsne_class =  clust10_tree, 
               logPLT = log(PLT),
            rat_BM_HB = BM_BLAST/HB,
            rat_WBC_ANC = WBC/ANC,
            rat_WBC_HB = WBC/HB,
            rat_WBC_PLT = WBC/logPLT) 

clinical_tree_sub <- filter(clinical_tree, meta %in% meta50) %>%
     filter(!meta %in% c(51)) %>% ## Remove a cluster with bad classification
     select(-tsne_class)

library(rpart)

tree_mod <- rpart(meta ~ . , data = clinical_tree_sub, method = "class", 
 control = rpart.control(cp = 0.001))

clinical_tree_sub$tree_pred <- predict(tree_mod, clinical_tree_sub, type = "class")
table(Pred = clinical_tree_sub$tree_pred, meta = clinical_tree_sub$meta)
clinical_tree$tree_pred <- predict(tree_mod, clinical_tree, type = "class")

table(Pred = clinical_tree$tree_pred, tsne = clinical_tree$tsne_class)
table(Pred = clinical_tree$tree_pred, meta = clinical_tree$meta)

umap_pred <- clinical_tree$tree_pred
save(tree_mod, umap_pred, file = "results/clustering/umap_meta_cluster.Rdata")

table(cutree(meta_umap, h = 5), clust10_tree)
### n_neighbors = 15
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res <- umap(clust_mat)

png("figures/umap_2comp_15neigh.png")
plot(umap_res$layout)
dev.off()

umap_dist15 <- dist(umap_res$layout)
hclust15 <- hclust(umap_dist15)

sil_scores_tree15 <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist15, clusth = hclust15)
which.max(sil_scores_tree15)  + 1
png("figures/umap15_silhouette_tree.png")
plot(2:40, sil_scores_tree15, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

### n_neighbors = 50
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res_50n <- umap(clust_mat, n_neighbors = 50)

png("figures/umap_2comp_50neigh.png")
plot(umap_res_50n$layout)
dev.off()


umap_dist50 <- dist(umap_res_50n$layout)
hclust50 <- hclust(umap_dist50)

sil_scores_tree50 <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist50, clusth = hclust50)
which.max(sil_scores_tree50)  + 1
png("figures/umap50_silhouette_tree.png")
plot(2:40, sil_scores_tree50, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust50 <- cutree(hclust50, k = 17)

table(clust50, clust10_tree)



### n_neighbors = 50
### n_components = 2
### min_dist = 0.01
set.seed(27)
umap_res_50n_01d <- umap(clust_mat, n_neighbors = 50, min_dist = 0.01)

png("figures/umap_2comp_50neigh_01d.png")
plot(umap_res_50n_01d$layout)
dev.off()


umap_dist50_01d <- dist(umap_res_50n_01d$layout)
hclust50_01d <- hclust(umap_dist50_01d)

sil_scores_tree50_01d <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist50_01d, clusth = hclust50_01d)
which.max(sil_scores_tree50_01d)  + 1
png("figures/umap50_01d_silhouette_tree.png")
plot(2:40, sil_scores_tree50_01d, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust50_01d <- cutree(hclust50, k = 18)



### n_neighbors = 50
### n_components = 6
### min_dist = 0.1
set.seed(27)
umap_res_50n_6c <- umap(clust_mat, n_neighbors = 50, n_components = 6)

png("figures/umap_6comp_50neigh.png")
plot(umap_res_50n_6c$layout)
dev.off()


umap_dist50_6c <- dist(umap_res_50n_6c$layout)
hclust50_6c <- hclust(umap_dist50_6c)

sil_scores_tree50_6c <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist50_6c, clusth = hclust50_6c)
png("figures/umap50_6c_silhouette_tree.png")
plot(2:40, sil_scores_tree50_6c, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust50_6c <- cutree(hclust50, k = 11)

### n_neighbors = 50
### n_components = 6
### min_dist = 0.01
set.seed(27)
umap_res_50n_6c_01d <- umap(clust_mat, n_neighbors = 50, min_dist = 0.01, n_components = 6)

png("figures/umap_6comp_50neigh_01d.png")
plot(umap_res_50n_6c_01d$layout)
dev.off()


umap_dist50_6c_01d <- dist(umap_res_50n_6c_01d$layout)
hclust50_6c_01d <- hclust(umap_dist50_6c_01d)

sil_scores_tree50_6c_01d <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist50_6c_01d, clusth = hclust50_01d)
which.max(sil_scores_tree50_01d)  + 1
png("figures/umap50_01d_silhouette_tree.png")
plot(2:40, sil_scores_tree50_01d, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust50_01d <- cutree(hclust50, k = 18)



### n_neighbors = 5
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res_5n <- umap(clust_mat, n_neighbors = 5)

png("figures/umap_2comp_5neigh.png")
plot(umap_res_5n$layout)
dev.off()

### n_neighbors = 10
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res_10n <- umap(clust_mat, n_neighbors = 10)

png("figures/umap_2comp_10neigh.png")
plot(umap_res_10n$layout)
dev.off()


### n_neighbors = 100
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res_100n <- umap(clust_mat, n_neighbors = 100)

png("figures/umap_2comp_100neigh.png")
plot(umap_res_100n$layout)
dev.off()


umap_dist100 <- dist(umap_res_100n$layout)
hclust100 <- hclust(umap_dist100)

sil_scores_tree100 <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist100, clusth = hclust100)
which.max(sil_scores_tree100)  + 1
png("figures/umap100_silhouette_tree.png")
plot(2:40, sil_scores_tree100, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust100 <- cutree(hclust100, k = 22)

table(clust100, clust10_tree)


### n_neighbors = 30
### n_components = 2
### min_dist = 0.1
set.seed(27)
umap_res_30n <- umap(clust_mat, n_neighbors = 30)

png("figures/umap_2comp_30neigh.png")
plot(umap_res_30n$layout)
dev.off()


umap_dist30 <- dist(umap_res_30n$layout)
hclust30 <- hclust(umap_dist30)

sil_scores_tree30 <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist30, clusth = hclust30)
which.max(sil_scores_tree30)  + 1
png("figures/umap30_silhouette_tree.png")
plot(2:40, sil_scores_tree30, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust30 <- cutree(hclust30, k = 17)

table(clust30, clust10_tree)





### n_neighbors = 30
### n_components = 2
### min_dist = 0.01
set.seed(27)
umap_res_30n_01d <- umap(clust_mat, n_neighbors = 30, min_dist = 0.01)

png("figures/umap_2comp_30neigh_01dist.png")
plot(umap_res_30n_01d$layout)
dev.off()


umap_dist30_01d <- dist(umap_res_30n_01d$layout)
hclust30_01d <- hclust(umap_dist30_01d)

sil_scores_tree30_01d <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist30_01d, clusth = hclust30_01d)
which.max(sil_scores_tree30_01d)  + 1
png("figures/umap30_1d_silhouette_tree.png")
plot(2:40, sil_scores_tree30_01d, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust30_01d <- cutree(hclust30_01d, k = 18)

table(clust30_01d, clust10_tree)


### n_neighbors = 30
### n_components = 2
### min_dist = 0.005
set.seed(27)
umap_res_30n_005d <- umap(clust_mat, n_neighbors = 30, min_dist = 0.005)

png("figures/umap_2comp_30neigh_005dist.png")
plot(umap_res_30n_005d$layout)
dev.off()


umap_dist30_005d <- dist(umap_res_30n_005d$layout)
hclust30_005d <- hclust(umap_dist30_005d)

sil_scores_tree30_005d <- sapply(2:40, silhouette_score_tree, data_dist = umap_dist30_005d, clusth = hclust30_005d)
which.max(sil_scores_tree30_005d)  + 1
png("figures/umap30_005d_silhouette_tree.png")
plot(2:40, sil_scores_tree30_005d, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

clust30_005d <- cutree(hclust30_005d, k = 12)

table(clust30_005d, clust10_tree)
