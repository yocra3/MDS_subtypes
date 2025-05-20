#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define clusters of patients using t-SNE
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


#' Define number of t-SNEs 
#' 
#' Run a PCA first to define N PCs
pca_mat <- prcomp(clust_mat)
summary(pca_mat)
# Importance of components:
#                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
# Standard deviation     1.5781 1.2195 1.0466 0.9867 0.92530 0.87173 0.81677
# Proportion of Variance 0.2614 0.1561 0.1150 0.1022 0.08989 0.07978 0.07004
# Cumulative Proportion  0.2614 0.4176 0.5326 0.6348 0.72466 0.80444 0.87448
#                            PC8     PC9    PC10    PC11    PC12   PC13    PC14
# Standard deviation     0.77487 0.48452 0.32430 0.26592 0.21931 0.1928 0.17372
# Proportion of Variance 0.06303 0.02465 0.01104 0.00742 0.00505 0.0039 0.00317
# Cumulative Proportion  0.93751 0.96216 0.97320 0.98062 0.98567 0.9896 0.99274
#                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
# Standard deviation     0.16491 0.09718 0.09591 0.09153 0.08968 0.07423 0.03664
# Proportion of Variance 0.00286 0.00099 0.00097 0.00088 0.00084 0.00058 0.00014
# Cumulative Proportion  0.99560 0.99659 0.99756 0.99844 0.99928 0.99986 1.00000

## Check 2, 6 and 8 PCs (41, 80 and 93%)

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

## Perpelixity 10 - 2 comps
set.seed(27)
tsne_res10 <- Rtsne(clust_mat, pca = FALSE,check_duplicates = FALSE, perplexity = 10) 

png("figures/tsne_clust10.png")
plot(tsne_res10$Y)
dev.off()

tsne_dist <- dist(tsne_res10$Y)
hclust10 <- hclust(tsne_dist)

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
sil_scores <- sapply(2:40, silhouette_score, data = tsne_res10$Y)
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


colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
             "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
             "#c5b0d5", "#c49c94", "#f7b6d2")


pdf("figures/tsne_clust10.pdf")
tibble(t1 = tsne_res10$Y[, 1],
    t2 = tsne_res10$Y[, 2],
    clust = factor(clust10_tree)) %>%
    ggplot(aes(x = t1, y = t2, col = clust)) +
        geom_point() +
        theme_bw() +
        scale_color_manual(values = colors) +
        ggtitle("Hclust")
tibble(t1 = tsne_res10$Y[, 1],
    t2 = tsne_res10$Y[, 2],
    clust = factor(clust10_k)) %>%
    ggplot(aes(x = t1, y = t2, col = clust)) +
        geom_point() +
        scale_color_manual(values = colors) +
        theme_bw() +
        ggtitle("Kmeans")
dev.off()

names(kar_events) <- kar_events
lapply(kar_events, function(x) table(clust_dataset[[x]], clust10_k))
lapply(kar_events, function(x) table(clust_dataset[[x]], clust10_tree))

table(clust10_tree, clinical_all$WHO_2016)
save(clust10_tree, tsne_res10, file = "results/tsne_perp10.Rdata")


## Perpelixity 10 - 3 comps
set.seed(27)
tsne_res10_3comp <- Rtsne(clust_mat, pca = FALSE, 
  check_duplicates = FALSE, perplexity = 10, dims = 3) 

png("figures/tsne_clust10_3comp.png")
plot(tsne_res10_3comp$Y)
dev.off()


tsne_dist_3comp <- dist(tsne_res10_3comp$Y)
hclust10_3comp <- hclust(dist(tsne_res10_3comp$Y))

sil_scores_tree_3comp <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist_3comp, clusth = hclust10_3comp)
png("figures/tsne10_3comp_silhouette_tree.png")
plot(2:40, sil_scores_tree_3comp, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree_3comp) + 1

clust10_3comp_tree <- cutree(hclust10_3comp, k = 15)
table(clust10_3comp_tree, clust10_tree)

png("figures/tsne_clust10_3comp_class.png")
tibble(t1 = tsne_res10$Y[, 1],
    t2 = tsne_res10$Y[, 2],
    clust = factor(clust10_3comp_tree)) %>%
    ggplot(aes(x = t1, y = t2, col = clust)) +
        geom_point() +
        scale_color_manual(values = colors) +
        theme_bw() 
dev.off()


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
clust30_tree <- cutree(hclust10_3comp, k = 16)

## Perplexity 30 - 3comp
set.seed(27)
tsne_res30_3comp <- Rtsne(clust_mat, pca = FALSE, 
  check_duplicates = FALSE, perplexity = 30, dims = 3) 

tsne_dist30_3comp <- dist(tsne_res30_3comp$Y)
hclust30_3comp <- hclust(tsne_dist30_3comp)

sil_scores_tree30_3comp <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist30_3comp, clusth = hclust30_3comp)
png("figures/tsne30_3comp_silhouette_tree.png")
plot(2:40, sil_scores_tree30_3comp, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree30_3comp) + 1

clust30_3comp_tree <- cutree(hclust30_3comp, k = 16)
table(clust30_3comp_tree, clust10_tree)
table(clust30_3comp_tree, clust10_3comp_tree)

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


## Perplexity 50 - 3comp
set.seed(27)
tsne_res50_3comp <- Rtsne(clust_mat, pca = FALSE, 
  check_duplicates = FALSE, perplexity = 50, dims = 3) 

tsne_dist50_3comp <- dist(tsne_res50_3comp$Y)
hclust50_3comp <- hclust(tsne_dist50_3comp)

sil_scores_tree50_3comp <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist50_3comp, clusth = hclust50_3comp)
png("figures/tsne50_3comp_silhouette_tree.png")
plot(2:40, sil_scores_tree50_3comp, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree50_3comp) + 1

clust50_3comp_tree <- cutree(hclust50_3comp, k = 15)
table(clust50_3comp_tree, clust10_tree)
table(clust50_3comp_tree, clust10_3comp_tree)

save(clust50_3comp_tree, clust10_3comp_tree, clust30_tree,  
     clust30_3comp_tree, file = "results/clustering/tsne_clusters.Rdata")

lapply(2:16, function(k){
    table(cutree(hclust10, k = k), clust10_tree)
})

 lapply(a, function(x) { sapply(seq_len(nrow(x)), function(i) {
 v <- x[i, ]
 v <- v[v != 0]
 paste(c(i, ":", names(v)), collapse = "_")
 })
})

## Perplexity 50 - 3comp
set.seed(27)
tsne_res50_3comp <- Rtsne(clust_mat, pca = FALSE, 
  check_duplicates = FALSE, perplexity = 50, dims = 3) 

tsne_dist50_3comp <- dist(tsne_res50_3comp$Y)
hclust50_3comp <- hclust(tsne_dist50_3comp)

sil_scores_tree50_3comp <- sapply(2:40, silhouette_score_tree, data_dist = tsne_dist50_3comp, clusth = hclust50_3comp)
png("figures/tsne50_3comp_silhouette_tree.png")
plot(2:40, sil_scores_tree50_3comp, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Silhouette Analysis")
dev.off()
which.max(sil_scores_tree50_3comp) + 1

clust50_3comp_tree <- cutree(hclust50_3comp, k = 16)
table(clust50_3comp_tree, clust10_tree)
table(clust50_3comp_tree, clust10_3comp_tree)
