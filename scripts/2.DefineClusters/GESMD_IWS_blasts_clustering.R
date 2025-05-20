#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide MDS blasts using PCA, CA and k-means
#' 
#' ---------------------------
#'
#' Notes:
#' Combine GESMD and IWS data. Compute PCA for clinical variables and CA for karyotype and mutation
#' data. Combine the components explaning >80% of variance in a new
#' matrix and run PCA again. Run k-means clustering on the final matrix.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(MASS)
library(cluster)
library(tidyverse)
library(rpart)
library(ca)
library(corrplot)


load("results/preprocess/clinical_preproc.Rdata")
load("data/GESMD/gesmd_data_all.Rdata")

mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Merge datasets
gesmd <- gesmd_data %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

com_variables <- intersect(colnames(gesmd), colnames(clinical))
## Select variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p", "complex",
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[!mutations %in% c("ID", "MLL_PTD", "ETNK1", "PPM1D", "BCORL1", "GNB1", "PRPF8")]

sel_vars <- c("ID", clin_vars, kar_events, mut_vars, "consensus", "SF3B1", "del5q", "WHO_2016")
sel_vars <- sel_vars[sel_vars %in% com_variables]

# Select data for clustering
comb_dataset <- rbind(select(gesmd, com_variables) %>% mutate(dataset = "GESMD"),
                       select(clinical, com_variables) %>% mutate(dataset = "IWS"))
clust_dataset <- comb_dataset %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    select(any_of(c("ID", clin_vars, kar_events, mut_vars, "dataset"))) %>%
    filter(complete.cases(.)) %>%
    filter(WBC < 20) %>% ## Remove outliers
    mutate(complex = ifelse(complex == "complex", 1, 0)) 

# Clinical PCA
clinvars_mat <- clust_dataset %>%
    select(BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
    as.matrix()

clin_pc <- prcomp(clinvars_mat, scale = TRUE)

df_clin <- as_tibble(clin_pc$x) %>%
    mutate(ID = clust_dataset$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID") 

png("figures/GESMD_IWS_clustering/clinical_PCs.png", width = 800, height = 800)
ggplot(df_clin, aes(x = PC1, y = PC2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()

# Karyotype CA
kar_mat <- clust_dataset[, kar_events[kar_events %in% colnames(clust_dataset)]] %>% as.matrix()
kar_mat  <- kar_mat [, colMeans(kar_mat ) > 0.01] ## Remove low frequency karyotypes
kar_mat <- (kar_mat - 1) * -1
kar_pca <- ca(kar_mat)

df_kar <- as_tibble(kar_pca$rowcoord) %>%
    mutate(ID = clust_dataset$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID")

png("figures/GESMD_IWS_clustering/karyotipe_PCs.png", width = 800, height = 800)
 df_kar %>%
    ggplot(aes(x = Dim1, y = Dim2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()

# Mutation PCA
mut_mat <- clust_dataset[, mut_vars[mut_vars %in% colnames(clust_dataset)]] %>% as.matrix()
mut_mat  <- mut_mat [, colMeans(mut_mat ) > 0.05] ## Remove low frequency mutations & TET2
mut_mat <- mut_mat[, colnames(mut_mat) != "TET2"]
mut_mat <- (mut_mat - 1) * -1
mut_pca <- ca(mut_mat)

df_mut <- as_tibble(mut_pca$rowcoord) %>%
    mutate(ID = clust_dataset$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID") 

png("figures/GESMD_IWS_clustering/mutation_PCs.png", width = 800, height = 800)
df_mut %>%
    ggplot(aes(x = Dim1, y = Dim2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()


# Combine components
# Select components explaining ~80% of variation
comb_pc_mat <- cbind(clin_pc$x[, 1:3], kar_pca$rowcoord[, 1:4], mut_pca$rowcoord[, 1:8])
comb_pc <- prcomp(comb_pc_mat)

df_comb <- as_tibble(comb_pc$x) %>%
    mutate(ID = clust_dataset$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID")

png("figures/GESMD_IWS_clustering/combined_PCs.png", width = 500, height = 500)
df_comb %>%
    ggplot(aes(x = PC1, y = PC2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()



png("figures/GESMD_IWS_clustering/corr_PCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat), method = 'number')
dev.off()

png("figures/GESMD_IWS_clustering/corr_oriPCs_combPCs.png", width = 800, height = 800)
corrplot(cor(comb_pc_mat, comb_pc$x[, 1:11]), method = 'number')
dev.off()


## K-means clustering
### Define number of clusters based on silhouette score
### Select 11 components (81.5% of variance)
set.seed(27)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

comb_pc_sel <- comb_pc$x[, 1:11]
dist_comb_pc <- dist(comb_pc_sel, method = "euclidean")
sil_scores <- sapply(2:20, silhouette_score,
      mat = comb_pc_sel, data_dist = dist_comb_pc)

png("figures/GESMD_IWS_clustering/comb_kmeans_silhouette.png")
plot(2:20, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 10, 6
clust_comb <- kmeans(comb_pc_sel, centers = 10, nstart = 1000)$cluster
table(clust_comb, clust_dataset$dataset)

lapply(c(colnames(kar_mat), colnames(mut_mat)), function(var){
    print(var)
    print(table(clust_dataset[[var]], clust_comb))
}) %>% invisible()


label_class <- 
    ifelse(clust_dataset$delY == 1, "Y-",
        ifelse(clust_dataset$del7 == 1, "7-",
            ifelse(clust_dataset$del20q == 1, "del20q",
                ifelse(clust_dataset$plus8 == 1, "8+", 
                    ifelse(clust_dataset$DNMT3A == 1, "DNMT3A",
                        ifelse(clust_dataset$STAG2 == 1, "STAG2",
                            ifelse(clust_dataset$TET2bi == 1, "TET2 bi-allelic",
                                ifelse(clust_dataset$TET2other == 1, "TET2 monoallelic",
                                    ifelse(clust_dataset$WBC > 6, "Midly Leukopenic",                 
                                           "Highly Leukopenic")))))))))
label_class <- factor(label_class, levels = c("8+", "STAG2", "del20q", "DNMT3A", "Midly Leukopenic", "TET2 monoallelic", "Y-", "TET2 bi-allelic", "Highly Leukopenic", "7-"))                                           
table(label_class, clust_comb)

tree_dataset <- select(clust_dataset, 
    all_of(c(colnames(clinvars_mat), colnames(kar_mat), colnames(mut_mat)))) %>%
    mutate(clust = label_class)

tree_mod <- rpart(clust ~ . , data = tree_dataset, method = "class")
predictions <- predict(tree_mod, tree_dataset, type = "class")
table(predictions, clust_comb)
table(predictions, label_class)

png("figures/GESMD_IWS_clustering/combined_classification_tree_ori.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()

## Re-assing samples to clusters
round <- 1
tree_class <- clust_comb
while (sum(tree_class != predictions) > 0) {
    tree_class <- predictions
    tree_dataset$clust <- predictions
    tree_mod <- rpart(clust ~ . , data = tree_dataset, method = "class")
    predictions <- predict(tree_mod, tree_dataset, type = "class")
    round <- round + 1
    print(round)
}
table(predictions, label_class)


png("figures/GESMD_IWS_clustering/combined_classification_tree.png", width = 1000, height = 650)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()

lapply(c(colnames(kar_mat), colnames(mut_mat)), function(var){
    print(var)
    print(table(clust_dataset[[var]], predictions))
}) %>% invisible()


# ## Relabel clusters and simplify classification
# tree_data_slim <- tree_dataset %>% 
#     mutate(clust = factor(labels[as.numeric(clust)], 
#            levels = c("Highly Leukopenic", "Mildly Leukopenic", "Y-", "7-", "del20q", "STAG2/8+"))) %>%
#     select("WBC", "STAG2", "del20q", "plus8", "del7", "delY", clust)

# tree_labels <- rpart(clust ~ . , data = tree_data_slim, method = "class")


# png("figures/GESMD_IWS_clustering/blasts_classification_tree_labels.png", width = 1000, height = 650)
# plot(tree_labels)
# text(tree_labels, use.n = TRUE)
# dev.off()
# table(new = predict(tree_labels, tree_data_slim, type = "class"), ori = tree_data_slim$clust)

df_clin$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "DNMT3A", "STAG2"))
df_kar$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "DNMT3A", "STAG2"))
df_mut$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "DNMT3A", "STAG2"))
df_comb$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "DNMT3A", "STAG2"))

save(sil_scores, clin_pc, df_clin, kar_pca, df_kar,  mut_pca, df_mut, comb_pc, df_comb, 
    file = "results/clustering/combined_pc_clustering_raw.Rdata")

#clust_dataset$clust_tree <-  predict(tree_mod, tree_dataset, type = "class")
clust_dataset$clust_manual <- label_class
clust_dataset$clust_ori <- clust_comb

combined_full <- left_join(clust_dataset, 
    select(clinical, hma, transplant, consensus, "ID", "IPSSM", "IPSSM_SCORE", "IPSSR", "IPSSR_SCORE", "AGE", "SEX", complex, starts_with("OS"), starts_with("AMLt")), by = "ID") 
save(combined_full, file = "results/clustering/combined_pc_clustering.Rdata")



# ## Re-cluster remaining samples
# clust_dataset$clust <-  predict(tree_labels, clust_dataset, type = "class")

# clust_dataset_mini <- clust_dataset %>%
#     filter(clust %in% c("Highly Leukopenic")) %>%
#     select(-STAG2, -del20q, -plus8, -del7, -delY) 

# # Clinical PCA
# clinvars_mat_mini <- clust_dataset_mini %>%
#     select(BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
#     mutate(WBC = ifelse(WBC > 13, 13, WBC),
#            ANC = ifelse(ANC > 10, 10, ANC),
#            MONOCYTES = ifelse(MONOCYTES > 3, 3, MONOCYTES)) %>%
#     as.matrix()

# clin_pc_mini <- prcomp(clinvars_mat_mini, scale = TRUE)

# png("figures/GESMD_IWS_clustering/clinical_PCs_noMol.png", width = 800, height = 800)
# as_tibble(clin_pc_mini$x) %>%
#     mutate(ID = clust_dataset_mini$ID) %>%
#     left_join(comb_dataset %>% select(ID, dataset), by = "ID") %>%
#     ggplot(aes(x = PC1, y = PC2, color = dataset)) +
#     geom_point() +
#     theme_bw()
# dev.off()


# # Mutation PCA
# mut_mat_mini <- clust_dataset_mini[, mut_vars[mut_vars %in% colnames(clust_dataset_mini)]] %>% as.matrix()
# mut_mat_mini  <- mut_mat_mini [, colMeans(mut_mat_mini ) > 0.05] ## Remove low frequency mutations & TET2
# mut_mat_mini <- mut_mat_mini[, colnames(mut_mat_mini) != "TET2"]
# mut_mat_mini <- (mut_mat_mini - 1) * -1
# mut_pca_mini <- ca(mut_mat_mini)

# # Combine components
# # Select components explaining ~80% of variation
# comb_pc_mat_mini <- cbind(clin_pc_mini$x[, 1:4], mut_pca_mini$rowcoord[, 1:6])
# comb_pc_mini <- prcomp(comb_pc_mat_mini)


# comb_pc_mini_sel <- comb_pc_mini$x[, 1:7]
# dist_comb_pc_mini <- dist(comb_pc_mini_sel, method = "euclidean")
# set.seed(27)
# sil_scores_mini <- sapply(2:20, silhouette_score,
#       mat = comb_pc_mini_sel, data_dist = dist_comb_pc_mini)

# png("figures/GESMD_IWS_clustering/comb_kmeans_silhouette_noMol.png")
# plot(2:20, sil_scores_mini, type="b", pch = 19, frame = FALSE,
#      xlab = "Number of clusters K",
#      ylab = "Silhouette Score",
#      main = "Combined clustering")
# dev.off()

# ## Best clusters:9,17
# clust_comb_mini <- kmeans(comb_pc_mini_sel, centers = 9, nstart = 1000)$cluster
# table(clust_comb_mini, clust_dataset_mini$dataset)
# lapply(colnames(mut_mat_mini), function(var){
#     print(var)
#     print(table(clust_dataset_mini[[var]], clust_comb_mini))
# }) %>% invisible()

# clust_comb_mini13 <- kmeans(comb_pc_mini_sel, centers = 13, nstart = 1000)$cluster
# table(clust_comb_mini13, clust_dataset_mini$dataset)
# lapply(colnames(mut_mat_mini), function(var){
#     print(var)
#     print(table(clust_dataset_mini[[var]], clust_comb_mini13))
# }) %>% invisible()





# ## Explore the effect of K on clustering
# ks <- 2:11
# clinical_low$complex <- factor(clinical_low$complex, levels = c("non-complex", "complex"))

# lapply(ks, function(k){
#     clust_comb <- kmeans(comb_pc$x[, 1:8], centers = k, nstart = 1000)$cluster
#     print(table(clust_comb, clinical_low$clust))
#     kars <- c(kar_events, "complex")
#     kar_tabs <- lapply(kars, function(x) table(clinical_low[[x]], clust_comb))
#     a <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
#     colnames(a) <- kars
#     print(a)

# })
