#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide of MDS clusters using t-SNE and mutation data
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
library(Rumap)
library(cluster)
library(umap)


load("results/preprocess/clinical_preproc.Rdata")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select variables for clustering
#' - Demographic: age, sex
#' - Clinical: Blood cell proportions
#' - Karyotype events: number of aberrations 
demo_vars <- c("AGE", "Sex")
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p",
                "del20q", "del7q", "plus8", "plus19", "del7")
var_sum <- colnames(clinical)[grep("N_", colnames(clinical))]
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

# Convert variables for clustering
#' Convert sex to dummy
#' Convert PLT to log
#' Select complete cases
#' Scale continous variables
clust_dataset <- mutate(clinical, Sex = ifelse(SEX == "F", 1, 0),
        PLT = log(PLT)) %>% 
    select(c(ID, Sex, AGE, clin_vars, kar_events, var_sum, mut_vars, consensus)) %>%
    filter(complete.cases(.))

cont_mat <- clust_dataset[, c("AGE", clin_vars, var_sum)] %>%
    data.matrix()
cont_mat_scale <- scale(cont_mat)

clust_mat <- cbind(cont_mat_scale,
    clust_dataset[, c("Sex", kar_events, mut_vars)] %>% data.matrix())
rownames(clust_mat) <- clust_dataset$ID

#' Run t-SNE per subtype
mdstypes <- unique(clust_dataset$consensus)
names(mdstypes) <- mdstypes

set.seed(27)
umap_types <- lapply(mdstypes, function(mds){

    mds_ids <- subset(clust_dataset, consensus == mds)$ID
    mat <- clust_mat[mds_ids, ]
    umap_res <- umap(mat, n_neighbors = 30, n_components = 2, min_dist = 0.3)
})

## Compute silhouettes for each subtype
silhouette_score_tree <- function(data_dist, clusth, k) {
  km <- cutree(clusth, k = k)
  ss <- silhouette(km, data_dist)
  mean(ss[, 3])
}

silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

compute_silhouette <- function(umap_res){

    umap_dist <- dist(umap_res$layout)
    hclust <- hclust(umap_dist)
    sil_scores_tree <- sapply(2:40, silhouette_score_tree,
     data_dist = umap_dist, clusth = hclust)
}
silhouettes <- lapply(umap_types, compute_silhouette)

compute_silhouette2 <- function(umap_res){
   
   umap_dist <- dist(umap_res$layout)
   sil_scores <- sapply(2:40, silhouette_score,
     mat = umap_res$layout, data_dist = umap_dist)
}

silhouettes2 <- lapply(umap_types, compute_silhouette2)


pdf("figures/umap_subtypes/umap_10_silhouette.pdf")
lapply(mdstypes, function(mds){
    plot(2:40, silhouettes2[[mds]], type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = mds)
} )
dev.off()

sel_clusters <- c(6, 4, 2, 3, 4, 3)
names(sel_clusters) <- mdstypes

clusters <- lapply(mdstypes, function(mds){
    k <- sel_clusters[mds]
    umap_res <- umap_types[[mds]]
    # umap_dist <- dist(umap_res$layout)
    # hclust <- hclust(umap_dist)

    # clust_class <- cutree(hclust, k = k)
    km <- kmeans(umap_res$layout, centers = k, nstart = 1000)
    clust_class <- km$clust
    paste(mds, clust_class, sep = "_")
})

pdf("figures/umap_subtypes/umap_10.pdf")
lapply(mdstypes, function(mds){

    df <- as_tibble(umap_types[[mds]]$layout) %>%
        mutate(col = clusters[[mds]])
    ggplot(df, aes(x = V1, y = V2, col = col)) +
        geom_point() +
        theme_bw() +
        ggtitle(mds)
} )
dev.off()


## Generate merged dataset with clustering data
clust_comb <- lapply(mdstypes, function(mds){

    filter(clust_dataset, consensus == mds) %>%
        mutate(cluster = clusters[[mds]])
}) %>%
    Reduce(f = rbind)

clinical_clust <- left_join(clinical, 
    select(clust_comb, ID, cluster), by = "ID")
save(clinical_clust, file = "results/clustering/umap_mutations_clustering.Rdata")
save(umap_types, file = "results/clustering/umap_mutations.Rdata")

