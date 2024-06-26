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
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.2 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(missMDA)
library(Rtsne)

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
kar_events <- c("monoy", "del5q", "del20q", "triso8", "mono7")


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

set.seed(27)
tsne_res30 <- Rtsne(clust_mat, pca = FALSE,check_duplicates = FALSE, perplexity = 30) 

png("figures/tsne_clust30.png")
plot(tsne_res30$Y)
dev.off()


set.seed(27)
tsne_res50 <- Rtsne(clust_mat, pca = FALSE,check_duplicates = FALSE, perplexity = 50) 

png("figures/tsne_clust50.png")
plot(tsne_res50$Y)
dev.off()

set.seed(27)
tsne_res10 <- Rtsne(clust_mat, pca = FALSE,check_duplicates = FALSE, perplexity = 10) 

png("figures/tsne_clust10.png")
plot(tsne_res10$Y)
dev.off()

hclust10 <- hclust(dist(tsne_res10$Y))
clust10 <- cutree(hclust10, k = 5)

png("figures/tsne_clust10.png")
plot(tsne_res10$Y, col = clust10)
dev.off()

png("figures/tsne_clust50.png")
plot(tsne_res50$Y, col = clust10)
dev.off()

# Elbow Method
wcss <- function(data, k) {
  kmeans(data, centers = k, nstart = 25)$tot.withinss
}
k_values <- 1:30
wcss_values <- sapply(k_values, wcss, data = tsne_res10$Y)

png("figures/tsne10_elbow.png")
plot(k_values, wcss_values, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares",
     main = "Elbow Method")
dev.off()

# Silhouette Analysis
silhouette_score <- function(data, k) {
  km <- kmeans(data, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(data))
  mean(ss[, 3])
}
sil_scores <- sapply(2:40, silhouette_score, data = dist(tsne_res10$Y))
png("figures/tsne10_silhouette.png")
plot(2:40, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()

which.max(sil_scores) + 1


silhouette_score_tree <- function(data_dist, clusth, k) {
  km <- cutree(clusth, k = k)
  ss <- silhouette(km, data_dist)
  mean(ss[, 3])
}

tsne_dist <- dist(tsne_res10$Y)
sil_scores_tree <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist, clusth = hclust10)
png("figures/tsne10_silhouette_tree.png")
plot(2:40, sil_scores_tree, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree) + 1

# Gap Statistic
gap_stat <- clusGap(tsne_res10$Y, FUN = kmeans, nstart = 25, K.max = 40, B = 50)
fviz_gap_stat(gap_stat)

clust10_tree <- cutree(hclust10, k = which.max(sil_scores_tree) + 1)
clust10_k <-  kmeans(tsne_res10$Y, centers = which.max(sil_scores) + 1, nstart = 100)$cluster

pdf("figures/tsne_clust10.pdf")
tibble(t1 = tsne_res10$Y[, 1],
    t2 = tsne_res10$Y[, 2],
    clust = factor(clust10_tree)) %>%
    ggplot(aes(x = t1, y = t2, col = clust)) +
        geom_point() +
        theme_bw() +
        ggtitle("Hclust")
tibble(t1 = tsne_res10$Y[, 1],
    t2 = tsne_res10$Y[, 2],
    clust = factor(clust10_k)) %>%
    ggplot(aes(x = t1, y = t2, col = clust)) +
        geom_point() +
        theme_bw() +
        ggtitle("Kmeans")
dev.off()

names(kar_events) <- kar_events
lapply(kar_events, function(x) table(clust_dataset[[x]], clust10_k))
lapply(kar_events, function(x) table(clust_dataset[[x]], clust10_tree))

table(clust10_tree, clinical_all$WHO_2016)
save(clust10_tree, tsne_res10, file = "results/tsne_perp10.Rdata")

## Perplexity 30
tsne_dist30 <- dist(tsne_res30$Y)
hclust30 <- hclust(dist(tsne_res30$Y))

sil_scores_tree30 <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist30, clusth = hclust30)
png("figures/tsne30_silhouette_tree.png")
plot(2:40, sil_scores_tree30, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree30) + 1


## Perplexity 50
tsne_dist50 <- dist(tsne_res50$Y)
hclust50 <- hclust(dist(tsne_res50$Y))

sil_scores_tree50 <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist50, clusth = hclust50)
png("figures/tsne50_silhouette_tree.png")
plot(2:40, sil_scores_tree50, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree50) + 1

lapply(2:16, function(k){
    table(cutree(hclust10, k = k), clust10_tree)
})

 lapply(a, function(x) { sapply(seq_len(nrow(x)), function(i) {
 v <- x[i, ]
 v <- v[v != 0]
 paste(c(i, ":", names(v)), collapse = "_")
 })
})
