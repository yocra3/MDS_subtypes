#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define clusters of patients by grouping different t-SNE clusters
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
library(cluster)

load( "results/tsne_perp10.Rdata")
load("results/clustering/tsne_clusters.Rdata")


meta_mat <- Reduce(cbind, 
    list(clust10_3comp_tree, clust10_tree, 
        clust30_3comp_tree, clust30_tree, clust50_3comp_tree))

dist_mat <- apply(meta_mat, 1, function(x) {
    apply(meta_mat, 1, function(y) sum(x != y))
})

meta_hclust <- hclust(as.dist(dist_mat))
table(cutree(meta_hclust, h = 2), clust10_tree)

meta_class <- cutree(meta_hclust, h = 1)
tb <- tibble(meta = meta_class, tsne10 = clust10_tree)
n_meta <- group_by(tb, meta) %>%
     summarize(n = n()) %>%
     filter(n >= 90)

tb_filt <- mutate(tb, meta2 = ifelse(!meta %in% n_meta$meta, NA, meta))
table(tb_filt$meta2, tb_filt$tsne10, useNA = "ifany")

tb <- mutate(tb, meta3 = ifelse(tsne10 == 3 & meta == 5, 
    paste0(meta, "m"), tsne10))

clust10_meta <- tb$meta3
save(clust10_meta, file = "results/clustering/tsne10b_clusters.Rdata") 