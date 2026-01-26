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
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(FactoMineR)
library(survival)
library(survminer)
library(cluster)
library(rpart)

load("results/preprocess/clinical_preproc.Rdata")
load("data/GESMD/gesmd_data_1125.Rdata")

mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Merge datasets
gesmd <- gesmd_data_1125 %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

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


IWS_dataset <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>% ## Filtro IWS
     filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    select(all_of(c(mut_vars, kar_events, clin_vars, "ID"))) %>%
    filter(complete.cases(.)) 


## Group genes by function
gene_functions <- read_csv("data/gene_functions2.csv")

gene_cats <- unique(gene_functions$Specific_Function)

## Group low frequency mutations (<4%) 
mut_freq <- colMeans(IWS_dataset[, mut_vars])
gene_functions$Frequency <- mut_freq[gene_functions$Gene]
low_freq_vars <- names(mut_freq[mut_freq < 0.04])

for (cat in gene_cats) {
    genes_in_cat <- gene_functions %>%
        filter(Specific_Function == cat & Frequency < 0.04) %>%
        pull(Gene) %>%
        intersect(colnames(IWS_dataset))
    if (length(genes_in_cat) > 1) {
        IWS_dataset[[cat]] <- ifelse(rowSums(IWS_dataset[, genes_in_cat], na.rm = TRUE) > 0, 1, 0)
    }
}

## Remove low frequency karyotype events 
kar_freq <- colMeans(IWS_dataset[, kar_events])
low_freq_vars <- c(low_freq_vars, names(kar_freq[kar_freq < 0.01]))

## Remove two categories with low frequency after grouping
low_freq_vars <- c(low_freq_vars, "RNA Helicases", "Receptor Tyrosine Kinases")

IWS_cluster <- IWS_dataset %>%
    select(-all_of(c(low_freq_vars, "ID", "TP53", "TET2"))) %>%
    mutate(across(-all_of(clin_vars), as_factor))



avail_vars_list <- lapply(list(clin_vars, kar_events, c(mut_vars, gene_cats)),
    function(x) intersect(x, colnames(IWS_cluster))) 
names(avail_vars_list) <- c("clin_vars", "kar_events", "genes_mut")
IWS_order <- unlist(avail_vars_list)


## Run Multiple Factor Analysis (MFA) with FactoMiner
IWS_mfa <- MFA(IWS_cluster[,IWS_order], 
    group = c(length(clin_vars), length(IWS_order) - length(clin_vars)),
    type = c("s", "n"), graph = FALSE, ncp = 7)



png("figures/GESMD_IWS_clustering/IWS_MFA_eigen.png")
plot(IWS_mfa$global.pca$eig[, 2])
abline(h = 1/length(IWS_order)*100, col = "red", lty = 2)
dev.off()

## Option 1: 7 dimensions (elbow method)
## K-means clustering
### Define number of clusters based on silhouette score
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

IWS_MFA_sel1 <- IWS_mfa$ind$coord[, 1:7]
dist_comb_IWS_pc1 <- dist(IWS_MFA_sel1, method = "euclidean")
set.seed(27)
sil_scores1 <- sapply(2:20, silhouette_score,
      mat = IWS_MFA_sel1, data_dist = dist_comb_IWS_pc1)
png("figures/GESMD_IWS_clustering/IWS_MFA1_kmeans_silhouette.png")
plot(2:20, sil_scores1, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 8
clust_comb1 <- kmeans(IWS_MFA_sel1, centers = 8, nstart = 1000)$cluster
tabs1 <- lapply(unlist(avail_vars_list[-1]), function(var){
    print(var)
    print(table(IWS_cluster[[var]], clust_comb1))
}) 
names(tabs1) <- unlist(avail_vars_list[-1])

# cl1
table(IWS_cluster[["PHF6"]] == 1 &  IWS_cluster[["Chromatin Histone Modifiers"]] == 1,  clust_comb1)
table(IWS_cluster[["PHF6"]] == 1 | IWS_dataset[["ARID2"]] == 1,  clust_comb1)
table(IWS_cluster[["U2AF1"]] == 1 ,  clust_comb1)
table(IWS_cluster[["PHF6"]] == 1 ,  clust_comb1)
table(IWS_cluster[["PHF6"]] == 1 | IWS_dataset[["U2AF1"]] == 1,  clust_comb1)

# cl2
table(IWS_cluster[["STAG2"]] == 1,  clust_comb1)

# cl3
table(IWS_cluster[["EZH2"]] == 1,  clust_comb1)


# cl5
table(IWS_cluster[["TET2bi"]] == 1,  clust_comb1)

# cl6 
table(IWS_cluster[["Other Functions"]] == 1,  clust_comb1)
table(IWS_cluster[["RAS MAPK Pathway"]] == 1,  clust_comb1)

# cl7
table(IWS_cluster[["SETBP1"]] == 1 |  IWS_cluster[["del7"]] == 1,  clust_comb1)



table(IWS_cluster[["NPM1"]] == 1 | IWS_cluster[["WT1"]] == 1, clust_comb2)
table(IWS_cluster[["ASXL1"]] == 1 & IWS_cluster[["STAG2"]] == 1, clust_comb2)
table(IWS_cluster[["SETBP1"]] == 1 |  IWS_cluster[["ETV6"]] == 1,  clust_comb2)
table(IWS_cluster[["TET2bi"]] == 1 | IWS_cluster[["ZRSR2"]] == 1, clust_comb2)

table(IWS_cluster[["U2AF1"]] == 1 & IWS_cluster[["BCOR"]] == 1, clust_comb2)


silhouette_score_HCPC <- function(mfa, data_dist, k) {
  clust <-  HCPC(mfa, nb.clust = k, graph = FALSE)
  ss <- silhouette(as.integer(clust$data.clust$clust), data_dist)
  mean(ss[, 3])
}

set.seed(27)
sil_scores_hc <- sapply(2:20, silhouette_score_HCPC,
      mfa = IWS_mfa, data_dist = dist_comb_IWS_pc1)
png("figures/GESMD_IWS_clustering/IWS_MFA1_hcpc_silhouette.png")
plot(2:20, sil_scores_hc, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()


hc_clust <- HCPC(IWS_mfa, nb.clust = 8, graph = FALSE)
comb_clusthc <- hc_clust$data.clust$clust
table(comb_clusthc, clust_comb1)
#             clust_comb1
# comb_clusthc   1   2   3   4   5   6   7   8
#            1   0   0   0   0   0   0   0 197
#            2   0   0   0   0 146   0   0   0
#            3   1   0   0 482   0   0   0   0
#            4  78   1   0   0   0   0   0   0
#            5   0   0   0   0   0  56   0   0
#            6   0   0   0   0   0   0  72   0
#            7   0   4  68   0   0   0   0   0
#            8   0 130   0   1   0   0   0   0
## Mismo clustering que kmeans

## GESMD clustering
gesmd_dataset <- gesmd %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    filter(WBC < 20 & ANC < 40) 

## Select GESMD variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' - Mutations
kars_gesmd <-  colnames(gesmd_dataset)[42:51]
muts_gesmd <- colnames(gesmd_dataset)[53:265]

## Remove two mutations not present in any sample
muts_gesmd <- muts_gesmd[!muts_gesmd %in% c("TP53loh", "MLL_PTD")]
## Group low frequency mutations (<4%) 
mut_freq_gesmd <- colMeans(gesmd_dataset[, muts_gesmd], na.rm = TRUE)
## Remove mutations with zero frequency
gene_functions$Frequency_gesmd <- mut_freq_gesmd[gene_functions$Gene]
low_freq_vars_gesmd <- names(mut_freq_gesmd[mut_freq_gesmd < 0.04])

for (cat in gene_cats) {
    genes_in_cat <- gene_functions %>%
        filter(Specific_Function == cat & Frequency_gesmd < 0.04) %>%
        pull(Gene) %>%
        intersect(colnames(gesmd_dataset))
    if (length(genes_in_cat) > 1) {
        gesmd_dataset[[cat]] <- ifelse(rowSums(gesmd_dataset[, genes_in_cat], na.rm = TRUE) > 0, 1, 0)
    }
}

## Remove low frequency karyotype events 
kar_freq_gesmd <- colMeans(gesmd_dataset[, kars_gesmd], na.rm = TRUE)
low_freq_vars_gesmd <- c(low_freq_vars_gesmd, names(kar_freq_gesmd[kar_freq_gesmd < 0.01]))

low_freq_vars_gesmd <- c(low_freq_vars_gesmd, "del7_7q") ## Remove category repeated

## Remove two categories with low frequency after grouping
low_freq_vars_gesmd <- c(low_freq_vars_gesmd, "General Signaling")
low_freq_vars_gesmd[is.na(low_freq_vars_gesmd)] <- "del17_17p"
gesmd_dataset_filt <- gesmd_dataset %>%
    select(any_of(c(clin_vars, kars_gesmd, muts_gesmd, gene_cats, "ID"))) %>%
    select(-all_of(c(low_freq_vars_gesmd, "TP53mut", "TET2", "TP53"))) %>%
    filter(complete.cases(.))


avail_vars_list_gesmd <- lapply(list(clin_vars, kars_gesmd, c(muts_gesmd, gene_cats)),
    function(x) intersect(x, colnames(gesmd_dataset_filt))) 
names(avail_vars_list_gesmd) <- c("clin_vars", "kar_events", "genes_mut")
gesmd_order <- unlist(avail_vars_list_gesmd)

GESMD_cluster <- gesmd_dataset_filt %>%
    select(-c("ID")) %>%
    mutate(across(-all_of(clin_vars), as_factor))

## Run Multiple Factor Analysis (MFA) with FactoMiner
GESMD_mfa <- MFA(GESMD_cluster[,gesmd_order], 
    group = c(length(clin_vars), length(gesmd_order) - length(clin_vars)),
    type = c("s", "n"), graph = FALSE, ncp = 9)

png("figures/GESMD_IWS_clustering/GESMD_MFA_eigen.png")
plot(GESMD_mfa$global.pca$eig[, 2])
abline(h = 1/length(gesmd_order)*100, col = "red", lty = 2)
dev.off()


GESMD_MFA_sel <- GESMD_mfa$ind$coord[, 1:9]
dist_comb_GESMD_pc1 <- dist(GESMD_MFA_sel, method = "euclidean")
set.seed(27)
sil_scores_gesmd <- sapply(2:20, silhouette_score,
      mat = GESMD_MFA_sel, data_dist = dist_comb_GESMD_pc1)
png("figures/GESMD_IWS_clustering/GESMD_MFA_kmeans_silhouette.png")
plot(2:20, sil_scores_gesmd, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 10
clust_gesmd <- kmeans(GESMD_MFA_sel, centers = 10, nstart = 1000)$cluster
tabs_gesmd <- lapply(unlist(avail_vars_list_gesmd[-1]), function(var){
    print(var)
    print(table(GESMD_cluster[[var]], clust_gesmd))
}) 
names(tabs_gesmd) <- unlist(avail_vars_list_gesmd[-1])

## cl 1
table(GESMD_cluster[["EZH2"]] == 1,  clust_gesmd)

## cl 2
table(GESMD_cluster[["del20q"]] == 1,  clust_gesmd)
table(GESMD_cluster[["U2AF1"]] == 1 | GESMD_cluster[["del20q"]] == 1,  clust_gesmd)

## cl 3
table(GESMD_cluster[["TET2bi"]] == 1,  clust_gesmd)


## cl 4
table(GESMD_cluster[["del7"]] == 1,  clust_gesmd)

## cl 8
table(GESMD_cluster[["RAS MAPK Pathway"]] == 1,  clust_gesmd)

# cl 9 - STAG2
table(GESMD_cluster[["STAG2"]] == 1,  clust_gesmd)
table(GESMD_cluster[["SRSF2"]] == 1,  clust_gesmd)
table(GESMD_cluster[["SRSF2"]] == 1 | GESMD_cluster[["STAG2"]] == 1,  clust_gesmd)
table((GESMD_cluster[["SRSF2"]] == 1 & GESMD_cluster[["ASXL1"]] == 1) | GESMD_cluster[["STAG2"]] == 1,  clust_gesmd)

## cl 10
table(GESMD_cluster[["CBL"]] == 1,  clust_gesmd)

set.seed(27)
sil_scores_hc_gesmd <- sapply(2:20, silhouette_score_HCPC,
      mfa = GESMD_mfa, data_dist = dist_comb_GESMD_pc1)
png("figures/GESMD_IWS_clustering/GESMD_MFA1_hcpc_silhouette.png")
plot(2:20, sil_scores_hc_gesmd, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()


hc_clust_gesmd <- HCPC(GESMD_mfa, nb.clust = 10, graph = FALSE)
gesmd_clusthc <- hc_clust_gesmd$data.clust$clust
table(gesmd_clusthc, clust_gesmd)


#              clust_gesmd
# gesmd_clusthc  1  2  3  4  5  6  7  8  9 10
#            1   0  0  0  7  0  0  0  0  0  0
#            2   0  0  0  0  0  0  0  0 24  0
#            3   0  0  0  0  0  0 42  0  0  0
#            4   0  0  0  0  0  0  0 11  0  0
#            5   0  0  0  0  0  0  0  0  0 17
#            6  22  0  0  0  0  0  0  0  1  0
#            7   0  0  0  0  0 94  0  0  1  0
#            8   0 24  0  0  0  0  0  0  0  0
#            9   0  0 26  0  0  0  0  0  0  0
#            10  0  0  0  0 29  0  0  0  1  0


hc_clust_gesmd2 <- HCPC(GESMD_mfa, nb.clust = 12, graph = FALSE)
gesmd_clusthc2 <- hc_clust_gesmd2$data.clust$clust
table(gesmd_clusthc2, clust_gesmd)
    

IWS_cluster_tree <- IWS_cluster %>%
    mutate(clust = as.factor(clust_comb1))
IWS_tree <- rpart(clust ~ ., data = IWS_cluster_tree, method = "class")

png("figures/GESMD_IWS_clustering/IWS_tree.png", width = 1000, height = 700)
plot(IWS_tree)
text(IWS_tree, use.n = TRUE)
dev.off()

gesmd_cluster_tree <- GESMD_cluster %>%
    mutate(clust = as.factor(clust_gesmd))
GESMD_tree <- rpart(clust ~ ., data = gesmd_cluster_tree, method = "class")

png("figures/GESMD_IWS_clustering/GESMD_tree.png", width = 1000, height = 700)
plot(GESMD_tree)
text(GESMD_tree, use.n = TRUE)
dev.off()



classifySamples <- function(df){
    new_class <- case_when(
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    df$del7 == 1       ~ "del7",
    df$STAG2 == 1       ~ "STAG2",
    df$PHF6 == 1        ~ "PHF6",
    TRUE                ~ "Other" 
  )
}
iws_class <- classifySamples(IWS_cluster)
iws_tab <- table(iws_class, clust_comb1)
matrixStats::rowMaxs(iws_tab) / rowSums(iws_tab)

gesmd_class <- classifySamples(GESMD_cluster)
gesmd_tab <- table(gesmd_class, clust_gesmd)
matrixStats::rowMaxs(gesmd_tab) / rowSums(gesmd_tab)

data.frame(group = names(matrixStats::rowMaxs(gesmd_tab)),
            IWS = matrixStats::rowMaxs(iws_tab) / rowSums(iws_tab),
           GESMD = matrixStats::rowMaxs(gesmd_tab) / rowSums(gesmd_tab)) %>%
    arrange(desc(IWS))

matrixStats::colMaxs(iws_tab) / colSums(iws_tab)
matrixStats::colMaxs(gesmd_tab) / colSums(gesmd_tab)

table(iws_class, clust_comb1)
#                     clust_comb1
# iws_class           1   2   3   4   5   6   7   8
#   del7              1   0   0   3   1   0  46   0
#   EZH2              5   3  51   4   4   6   5  10
#   Other            31  34   7 461  29  43  20 181
#   PHF6             33   0   0   0   1   0   0   0
#   STAG2             3  87  10  13   0   4   1   2
#   TET2 bi-allelic   6  11   0   2 111   3   0   4

table(gesmd_class, clust_gesmd)
#                     clust_gesmd
# gesmd_class        1  2  3  4  5  6  7  8  9 10
#   del7             0  0  0  7  0  0  0  0  0  1
#   EZH2            18  0  1  0  1  2  1  2  3  0
#   Other            4 15  6  0 23 88 34  9  9 15
#   PHF6             0  9  0  0  0  1  0  0  0  0
#   STAG2            0  0  0  0  0  0  5  0 15  0
#   TET2 bi-allelic  0  0 19  0  5  3  2  0  0  1

save(sil_scores1, IWS_mfa, clust_comb1, IWS_dataset, IWS_cluster,
     sil_scores_gesmd, GESMD_mfa, clust_gesmd, gesmd_dataset_filt, GESMD_cluster,
    file = "results/clustering/MFA_results.Rdata")


### Define first round of clusters:
#' - EZH2
#' - PHF6
#' - TET2 bi-allelic
#' - STAG2
#' - del7

## Not working well
## Run cluster again without these mutations
IWS_dataset2 <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>% ## Filtro IWS
    filter(SF3B1 == 0 & del5q == 0 & EZH2 == 0 & PHF6 == 0 & TET2bi == 0 & STAG2 == 0 & del7 == 0) %>% ## Remove previous clusters
    select(all_of(c(mut_vars, kar_events, clin_vars, "ID"))) %>%
    filter(complete.cases(.)) 



## Group low frequency mutations (<4%) 
gene_functions2 <- read_csv("data/gene_functions2.csv")

mut_freq2 <- colMeans(IWS_dataset2[, mut_vars])
gene_functions2$Frequency <- mut_freq2[gene_functions2$Gene]
low_freq_vars2 <- names(mut_freq2[mut_freq2 < 0.04])

for (cat in gene_cats) {
    genes_in_cat <- gene_functions2 %>%
        filter(Specific_Function == cat & Frequency < 0.04) %>%
        pull(Gene) %>%
        intersect(colnames(IWS_dataset2))
    if (length(genes_in_cat) > 1) {
        IWS_dataset2[[cat]] <- ifelse(rowSums(IWS_dataset2[, genes_in_cat], na.rm = TRUE) > 0, 1, 0)
    }
}

## Remove low frequency karyotype events 
kar_freq2 <- colMeans(IWS_dataset2[, kar_events])
low_freq_vars2 <- c(low_freq_vars2, names(kar_freq2[kar_freq2 < 0.01]))
## Remove two categories with low frequency after grouping
low_freq_vars2 <- c(low_freq_vars2, "RNA Helicases", "Receptor Tyrosine Kinases")

IWS_cluster2 <- IWS_dataset2 %>%
    select(-all_of(c(low_freq_vars2, "ID", "TP53", "TET2"))) %>%
    mutate(across(-all_of(clin_vars), as_factor))

avail_vars_list2 <- lapply(list(clin_vars, kar_events, c(mut_vars, gene_cats)),
    function(x) intersect(x, colnames(IWS_cluster2))) 
names(avail_vars_list2) <- c("clin_vars", "kar_events", "genes_mut")
IWS_order2 <- unlist(avail_vars_list2)

## Run Multiple Factor Analysis (MFA) with FactoMiner
IWS_mfa2 <- MFA(IWS_cluster2[,IWS_order2], 
    group = c(length(clin_vars), length(IWS_order2) - length(clin_vars)),
    type = c("s", "n"), graph = FALSE, ncp = 7)



png("figures/GESMD_IWS_clustering/IWS_MFA2_eigen.png")
plot(IWS_mfa2$global.pca$eig[, 2])
abline(h = 1/length(IWS_order2)*100, col = "red", lty = 2)
dev.off()

## Option 1: 7 dimensions (elbow method)
## K-means clustering
### Define number of clusters based on silhouette score
IWS_MFA_sel2 <- IWS_mfa2$ind$coord[, 1:7]
dist_comb_IWS_pc2 <- dist(IWS_MFA_sel2, method = "euclidean")
set.seed(27)
sil_scores2 <- sapply(2:20, silhouette_score,
      mat = IWS_MFA_sel2, data_dist = dist_comb_IWS_pc2)
png("figures/GESMD_IWS_clustering/IWS_MFA2_kmeans_silhouette.png")
plot(2:20, sil_scores2, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 8
clust_comb2 <- kmeans(IWS_MFA_sel2, centers = 8, nstart = 1000)$cluster
tabs2 <- lapply(unlist(avail_vars_list2[-1]), function(var){
    print(var)
    print(table(IWS_cluster2[[var]], clust_comb2))
}) 
names(tabs2) <- unlist(avail_vars_list2[-1])

IWS_cluster_tree2 <- IWS_cluster2 %>%
    mutate(clust = as.factor(clust_comb2))
IWS_tree2 <- rpart(clust ~ ., data = IWS_cluster_tree2, method = "class")

png("figures/GESMD_IWS_clustering/IWS_tree2.png", width = 1000, height = 700)
plot(IWS_tree2)
text(IWS_tree2, use.n = TRUE)
dev.off()



## GESMD clustering
gesmd_dataset2 <- gesmd %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(SF3B1 == 0 & del5q == 0 & TET2bi == 0 & PHF6 == 0 & STAG2 == 0 & EZH2 == 0 & del7 == 0) %>% ## Remove SF3B1 and del5q
    filter(WBC < 20 & ANC < 40) 

## Select GESMD variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' - Mutations

## Remove two mutations not present in any sample
muts_gesmd2 <- muts_gesmd[!muts_gesmd %in% c("TP53loh", "MLL_PTD")]
## Group low frequency mutations (<4%) 
mut_freq_gesmd2 <- colMeans(gesmd_dataset2[, muts_gesmd2], na.rm = TRUE)
## Remove mutations with zero frequency
gene_functions2$Frequency_gesmd <- mut_freq_gesmd2[gene_functions2$Gene]
low_freq_vars_gesmd2 <- names(mut_freq_gesmd2[is.na(mut_freq_gesmd2) | mut_freq_gesmd2 < 0.04])

for (cat in gene_cats) {
    genes_in_cat <- gene_functions2 %>%
        filter(Specific_Function == cat & Frequency_gesmd < 0.04) %>%
        pull(Gene) %>%
        intersect(colnames(gesmd_dataset2))
    if (length(genes_in_cat) > 1) {
        gesmd_dataset2[[cat]] <- ifelse(rowSums(gesmd_dataset2[, genes_in_cat], na.rm = TRUE) > 0, 1, 0)
    }
}

## Remove low frequency karyotype events 
kar_freq_gesmd2 <- colMeans(gesmd_dataset2[, kars_gesmd], na.rm = TRUE)
low_freq_vars_gesmd2 <- c(low_freq_vars_gesmd2, names(kar_freq_gesmd2[kar_freq_gesmd2 < 0.01]))

low_freq_vars_gesmd2 <- c(low_freq_vars_gesmd2, "del7_7q") ## Remove category repeated

## Remove two categories with low frequency after grouping
low_freq_vars_gesmd2 <- c(low_freq_vars_gesmd2, "General Signaling", "Splicing Factors Core", "ATG2B")
low_freq_vars_gesmd2[is.na(low_freq_vars_gesmd2)] <- "del17_17p"
gesmd_dataset_filt2 <- gesmd_dataset2 %>%
    select(any_of(c(clin_vars, kars_gesmd, muts_gesmd, gene_cats, "ID"))) %>%
    select(-all_of(c(low_freq_vars_gesmd2, "TP53mut", "TET2", "TP53"))) %>%
    filter(complete.cases(.))


avail_vars_list_gesmd2 <- lapply(list(clin_vars, kars_gesmd, c(muts_gesmd, gene_cats)),
    function(x) intersect(x, colnames(gesmd_dataset_filt2))) 
names(avail_vars_list_gesmd2) <- c("clin_vars", "kar_events", "genes_mut")
gesmd_order2 <- unlist(avail_vars_list_gesmd2)

GESMD_cluster2 <- gesmd_dataset_filt2 %>%
    select(-c("ID")) %>%
    mutate(across(-all_of(clin_vars), as_factor))

## Run Multiple Factor Analysis (MFA) with FactoMiner
GESMD_mfa2 <- MFA(GESMD_cluster2[,gesmd_order2], 
    group = c(length(clin_vars), length(gesmd_order2) - length(clin_vars)),
    type = c("s", "n"), graph = FALSE, ncp = 15)

png("figures/GESMD_IWS_clustering/GESMD_MFA2_eigen.png")
plot(GESMD_mfa2$global.pca$eig[, 2])
abline(h = 1/length(gesmd_order2)*100, col = "red", lty = 2)
dev.off()


GESMD_MFA2_sel <- GESMD_mfa2$ind$coord[, 1:8]
dist_comb_GESMD_pc2 <- dist(GESMD_MFA2_sel, method = "euclidean")
set.seed(27)
sil_scores_gesmd2 <- sapply(2:20, silhouette_score,
      mat = GESMD_MFA2_sel, data_dist = dist_comb_GESMD_pc2)
png("figures/GESMD_IWS_clustering/GESMD_MFA2_kmeans_silhouette.png")
plot(2:20, sil_scores_gesmd2, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 13
clust_gesmd2 <- kmeans(GESMD_MFA2_sel, centers = 8, nstart = 1000)$cluster
tabs_gesmd2 <- lapply(unlist(avail_vars_list_gesmd2[-1]), function(var){
    print(var)
    print(table(GESMD_cluster2[[var]], clust_gesmd2))
}) 
names(tabs_gesmd2) <- unlist(avail_vars_list_gesmd2[-1])

gesmd_cluster_tree2 <- GESMD_cluster2 %>%
    mutate(clust = as.factor(clust_gesmd2))
GESMD_tree2 <- rpart(clust ~ ., data = gesmd_cluster_tree2, method = "class")

png("figures/GESMD_IWS_clustering/GESMD_tree2.png", width = 1000, height = 700)
plot(GESMD_tree2)
text(GESMD_tree2, use.n = TRUE)
dev.off()