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
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
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
load("data/GESMD/gesmd_data_1125.Rdata")

mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Merge datasets
gesmd <- gesmd_data_1125 %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

com_variables <- intersect(colnames(gesmd), colnames(clinical))
## Perform clustering independently
## Select IWS variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p", 
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[!mutations %in% c("ID")]

gene_functions <- read_csv("data/gene_functions.csv")
# comb_dataset <- rbind(select(gesmd, com_variables) %>% mutate(dataset = "GESMD"),
#                        select(clinical, com_variables) %>% mutate(dataset = "IWS"))
IWS_dataset <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>% ## Filtro IWS
     filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    filter(WBC < 20 & ANC < 40) %>% ## Remove outliers
    mutate(complex = ifelse(complex == "complex", 1, 0)) %>%
    select(all_of(c(mut_vars, kar_events, clin_vars, "ID"))) %>%
    filter(complete.cases(.)) 

gene_cats <- unique(gene_functions$Specific_Function)
for (cat in gene_cats) {
    genes_in_cat <- gene_functions %>%
        filter(Specific_Function == cat) %>%
        pull(Gene) %>%
        intersect(colnames(IWS_dataset))
    if (length(genes_in_cat) > 1) {
        IWS_dataset[[cat]] <- ifelse(rowSums(IWS_dataset[, genes_in_cat], na.rm = TRUE) > 0, 1, 0)
    }
}

# Clinical PCA
clinvars_IWS <- IWS_dataset %>%
    select(BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
    as.matrix()

clin_IWS_pc <- prcomp(clinvars_IWS, scale = TRUE)

df_IWS_clin <- as_tibble(clin_IWS_pc$x) %>%
    mutate(ID = IWS_dataset$ID) 

png("figures/GESMD_IWS_clustering/clinical_IWS_PCs.png", width = 800, height = 800)
ggplot(df_IWS_clin, aes(x = PC1, y = PC2)) +
    geom_point() +
    theme_bw()
dev.off()


# ## Remove low frequency karyotipic events (<5%) and mutations (<5%) in dataset & TET2
# kar_events <- kar_events[kar_events %in% colnames(clust_dataset)]
# kar_events_sel <- kar_events[colMeans(clust_dataset[clust_dataset$dataset == "IWS", kar_events]) > 0.01]
# kar_events_sel <- c(kar_events_sel, "complex")

# mut_vars <- mut_vars[mut_vars %in% colnames(clust_dataset)]
# mut_vars_sel <- mut_vars[colMeans(clust_dataset[clust_dataset$dataset == "IWS", mut_vars], na.rm = TRUE) > 0.01]

# clust_dataset2 <- clust_dataset %>%
#     select(any_of(c("ID", clin_vars, kar_events_sel, mut_vars_sel, "dataset"))) %>%
#     filter(complete.cases(.)) 


# Karyotype and mutation CA
# ## Remove low frequency karyotipic events and mutations (<1%) 
kar_mut_IWS_mat <- IWS_dataset[, c(kar_events, mut_vars, gene_cats)] %>% as.matrix()
kar_mut_IWS_mat <- kar_mut_IWS_mat[, !colnames(kar_mut_IWS_mat) %in% c("TET2", "TP53")]
kar_mut_IWS_mat <- kar_mut_IWS_mat[, colMeans(kar_mut_IWS_mat) > 0.01]
kar_mut_IWS_mat <- (kar_mut_IWS_mat - 1) * -1
kar_mut_IWS_pca <- ca(kar_mut_IWS_mat)

df_kar_mut_IWS <- as_tibble(kar_mut_IWS_pca$rowcoord) %>%
    mutate(ID = IWS_dataset$ID)

png("figures/GESMD_IWS_clustering/karyotype_mutation_IWS_PCs.png", width = 800, height = 800)
 df_kar_mut_IWS %>%
    ggplot(aes(x = Dim1, y = Dim2)) +
    geom_point() +
    theme_bw()
dev.off()

cumsum(kar_mut_IWS_pca$sv)/sum(kar_mut_IWS_pca$sv)
plot(kar_mut_IWS_pca$sv/sum(kar_mut_IWS_pca$sv))

comb_IWS_pc_mat <- cbind(clin_IWS_pc$x[, 1:4], kar_mut_IWS_pca$rowcoord[, 1:11])
comb_IWS_pc <- prcomp(comb_IWS_pc_mat)

df_IWS_comb <- as_tibble(comb_IWS_pc$x) %>%
    mutate(ID = IWS_dataset$ID)

png("figures/GESMD_IWS_clustering/combined_IWS_PCs.png", width = 500, height = 500)
df_IWS_comb %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point() +
    theme_bw()
dev.off()



png("figures/GESMD_IWS_clustering/IWS_corr_PCs.png", width = 800, height = 800)
corrplot(cor(comb_IWS_pc_mat), method = 'number')
dev.off()

png("figures/GESMD_IWS_clustering/IWS_corr_oriPCs_combPCs.png", width = 800, height = 800)
corrplot(cor(comb_IWS_pc_mat, comb_IWS_pc$x), method = 'number')
dev.off()


## K-means clustering
### Define number of clusters based on silhouette score
### Select 13 components (84% of variance)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}


comb_IWS_pc_sel <- comb_IWS_pc$x[, 1:15]
dist_comb_IWS_pc <- dist(comb_IWS_pc_sel, method = "euclidean")
set.seed(27)
sil_scores <- sapply(2:20, silhouette_score,
      mat = comb_IWS_pc_sel, data_dist = dist_comb_IWS_pc)
png("figures/GESMD_IWS_clustering/IWS_kmeans_silhouette.png")
plot(2:20, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 11, 14
clust_comb <- kmeans(comb_IWS_pc_sel, centers = 6, nstart = 1000)$cluster

lapply(colnames(kar_mut_IWS_mat), function(var){
    print(var)
    print(table(kar_mut_IWS_mat[, var], clust_comb))
}) %>% invisible()








# Select data for clustering
comb_dataset <- rbind(select(gesmd, com_variables) %>% mutate(dataset = "GESMD"),
                       select(clinical, com_variables) %>% mutate(dataset = "IWS"))
clust_dataset <- comb_dataset %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>% ## Filtro IWS
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    filter(WBC < 20 & ANC < 40) %>% ## Remove outliers
    mutate(complex = ifelse(complex == "complex", 1, 0)) 

## Remove low frequency karyotipic events (<5%) and mutations (<5%) in dataset & TET2
kar_events <- kar_events[kar_events %in% colnames(clust_dataset)]
kar_events_sel <- kar_events[colMeans(clust_dataset[clust_dataset$dataset == "IWS", kar_events]) > 0.01]
kar_events_sel <- c(kar_events_sel, "complex")

mut_vars <- mut_vars[mut_vars %in% colnames(clust_dataset)]
mut_vars_sel <- mut_vars[colMeans(clust_dataset[clust_dataset$dataset == "IWS", mut_vars], na.rm = TRUE) > 0.01]

clust_dataset2 <- clust_dataset %>%
    select(any_of(c("ID", clin_vars, kar_events_sel, mut_vars_sel, "dataset"))) %>%
    filter(complete.cases(.)) 

# Clinical PCA
clinvars_mat <- clust_dataset2 %>%
    select(BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT) %>%
    as.matrix()

clin_pc <- prcomp(clinvars_mat, scale = TRUE)

df_clin <- as_tibble(clin_pc$x) %>%
    mutate(ID = clust_dataset2$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID") 

png("figures/GESMD_IWS_clustering/clinical_PCs.png", width = 800, height = 800)
ggplot(df_clin, aes(x = PC1, y = PC2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()

# Karyotype and mutation CA
kar_mut_mat <- clust_dataset2[, c(kar_events_sel, mut_vars_sel)] %>% as.matrix()
kar_mut_mat <- (kar_mut_mat - 1) * -1
kar_mut_pca <- ca(kar_mut_mat)

df_kar_mut <- as_tibble(kar_mut_pca$rowcoord) %>%
    mutate(ID = clust_dataset2$ID) %>%
    left_join(comb_dataset %>% select(ID, dataset), by = "ID")

png("figures/GESMD_IWS_clustering/karyotype_mutation_PCs.png", width = 800, height = 800)
 df_kar_mut %>%
    ggplot(aes(x = Dim1, y = Dim2, color = dataset)) +
    geom_point() +
    theme_bw()
dev.off()

# # Mutation PCA
# mut_mat <- clust_dataset[, mut_vars[mut_vars %in% colnames(clust_dataset)]] %>% as.matrix()
# mut_mat  <- mut_mat [, colMeans(mut_mat ) > 0.05] ## Remove low frequency mutations & TET2
# mut_mat <- mut_mat[, colnames(mut_mat) != "TET2"]
# mut_mat <- (mut_mat - 1) * -1
# mut_pca <- ca(mut_mat)

# df_mut <- as_tibble(mut_pca$rowcoord) %>%
#     mutate(ID = clust_dataset$ID) %>%
#     left_join(comb_dataset %>% select(ID, dataset), by = "ID") 

# png("figures/GESMD_IWS_clustering/mutation_PCs.png", width = 800, height = 800)
# df_mut %>%
#     ggplot(aes(x = Dim1, y = Dim2, color = dataset)) +
#     geom_point() +
#     theme_bw()
# dev.off()


# Combine components
# Select components explaining >80% of variation
# cumsum(kar_pca$sv)/sum(kar_pca$sv)
# cumsum(mut_pca$sv)/sum(mut_pca$sv)
# comb_pc_mat <- cbind(clin_pc$x[, 1:4], kar_pca$rowcoord[, 1:5], mut_pca$rowcoord[, 1:8])
# comb_pc <- prcomp(comb_pc_mat)
cumsum(kar_mut_pca$sv)/sum(kar_mut_pca$sv)
comb_pc_mat <- cbind(clin_pc$x[, 1:4], kar_mut_pca$rowcoord[, 1:20])
comb_pc <- prcomp(comb_pc_mat)

df_comb <- as_tibble(comb_pc$x) %>%
    mutate(ID = clust_dataset2$ID) %>%
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
corrplot(cor(comb_pc_mat, comb_pc$x[, 1:18]), method = 'number')
dev.off()


## K-means clustering
### Define number of clusters based on silhouette score
### Select 13 components (84% of variance)
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}


comb_pc_sel <- comb_pc$x[, 1:18]
dist_comb_pc <- dist(comb_pc_sel, method = "euclidean")
set.seed(27)
sil_scores <- sapply(2:20, silhouette_score,
      mat = comb_pc_sel, data_dist = dist_comb_pc)

png("figures/GESMD_IWS_clustering/comb_kmeans_silhouette.png")
plot(2:20, sil_scores, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 11, 14
clust_comb <- kmeans(comb_pc_sel, centers = 16, nstart = 1000)$cluster
table(clust_comb, clust_dataset2$dataset)

lapply(colnames(kar_mut_mat), function(var){
    print(var)
    print(table(clust_dataset2[[var]], clust_comb))
}) %>% invisible()


label_class <- 
    ifelse(clust_dataset$delY == 1, "Y-",
        ifelse(clust_dataset$del7 == 1, "7-",
        ifelse(clust_dataset$complex == 1, "complex",
            ifelse(clust_dataset$del20q == 1, "del20q",
                ifelse(clust_dataset$plus8 == 1, "8+", 
                    ifelse(clust_dataset$del7q == 1, "del7q",
                            ifelse(clust_dataset$TET2bi == 1, "TET2 bi-allelic",
                                ifelse(clust_dataset$TET2other == 1, "TET2 monoallelic",
                                    ifelse(clust_dataset$STAG2 == 1, "STAG2",            
                                    ifelse(clust_dataset$WBC > 6, "Midly Leukopenic",                 
                                           "Highly Leukopenic"))))))))))
label_class <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "del7q", "8+", "Y-", "complex", "TET2 monoallelic", "STAG2", "del20q", "7-", "TET2 bi-allelic"))
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


df_clin$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))
df_kar$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))
df_mut$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))
df_comb$clust <- factor(label_class, levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2"))

save(sil_scores, clin_pc, df_clin, kar_pca, df_kar,  mut_pca, df_mut, comb_pc, df_comb, 
    file = "results/clustering/combined_pc_clustering_raw.Rdata")

clust_dataset$clust_tree <-  predict(tree_mod, tree_dataset, type = "class")
clust_dataset$clust_manual <- label_class
clust_dataset$clust_ori <- clust_comb

combined_full <- left_join(clust_dataset %>% select(-PLT), 
    select(comb_dataset, consensus, "ID", "IPSSM",  "IPSSR", "AGE", "SEX", complex, "PLT"), by = "ID") 
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
