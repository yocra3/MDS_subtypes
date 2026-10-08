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
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
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
gesmd <- gesmd_data_1125_filt %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

## Perform clustering independently
## Select IWS variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p", "complex",
                "del20q", "del7q", "plus8", "plus19", "del7", "del13", 
                "del9q", "del3q")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[!mutations %in% c("ID")]


IWS_dataset <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("MDS-LB", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other") & MDS_TYPE == "primary") %>% ## Filtro IWS
    select(any_of(c(mut_vars, kar_events, clin_vars, "ID"))) %>%
    filter(complete.cases(.)) %>%
    mutate(complex = ifelse(complex == "complex", 1, 0))


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
    type = c("s", "n"), graph = FALSE, ncp = 12)

png("figures/GESMD_IWS_clustering/clustering/IWS_MFA_eigen.png")
plot(IWS_mfa$global.pca$eig[, 2])
abline(h = 1/length(IWS_order)*100, col = "red", lty = 2)
dev.off()

## Option 1: 9 dimensions (elbow method)
## K-means clustering
### Define number of clusters based on silhouette score
silhouette_score <- function(mat, data_dist, k) {
  km <- kmeans(mat, centers = k, nstart = 1000)
  ss <- silhouette(km$cluster, data_dist)
  mean(ss[, 3])
}

IWS_MFA_sel1 <- IWS_mfa$ind$coord[, 1:8]
dist_comb_IWS_pc1 <- dist(IWS_MFA_sel1, method = "euclidean")
set.seed(27)
sil_scores1 <- sapply(2:20, silhouette_score,
      mat = IWS_MFA_sel1, data_dist = dist_comb_IWS_pc1)
png("figures/GESMD_IWS_clustering/clustering/IWS_MFA1_kmeans_silhouette.png")
plot(2:20, sil_scores1, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 11
set.seed(27)
clust_comb1 <- kmeans(IWS_MFA_sel1, centers = 9, nstart = 1000)$cluster
tabs1 <- lapply(unlist(avail_vars_list[-1]), function(var){
    print(var)
    print(table(IWS_cluster[[var]], clust_comb1))
}) 
names(tabs1) <- unlist(avail_vars_list[-1])



# cl3
table(IWS_cluster[["Other Functions"]] == 1,  clust_comb1)

# cl5
table(IWS_cluster[["complex"]] == 1,  clust_comb1)


# cl1
table(IWS_cluster[["del5q"]] == 1,  clust_comb1)
table(IWS_cluster[["SF3B1"]] == 1,  clust_comb1)

table(IWS_cluster[["del5q"]] == 1 | IWS_cluster[["SF3B1"]] == 1,  clust_comb1)

# cl7
table(IWS_cluster[["EZH2"]] == 1,  clust_comb1)

# cl7
table(IWS_cluster[["del7"]] == 1,  clust_comb1)
table(IWS_cluster[["SETBP1"]] == 1,  clust_comb1)
table(IWS_cluster[["del7"]] == 1 | IWS_cluster[["SETBP1"]] == 1,  clust_comb1)


## cl4
table(IWS_cluster[["PHF6"]] == 1,  clust_comb1)
table(IWS_cluster[["RUNX1"]] == 1,  clust_comb1)

# cl2
table(IWS_cluster[["STAG2"]] == 1,  clust_comb1)
table(IWS_cluster[["SRSF2"]] == 1,  clust_comb1)
table(IWS_cluster[["IDH2"]] == 1,  clust_comb1)
table(IWS_cluster[["RUNX1"]] == 1,  clust_comb1)

table(IWS_cluster[["STAG2"]] == 1  &  IWS_cluster[["EZH2"]] == 0,  clust_comb1)


table((IWS_cluster[["STAG2"]] == 1 &  IWS_cluster[["SRSF2"]] == 1 )| IWS_cluster[["IDH2"]] == 1,  clust_comb1)
table(IWS_cluster[["STAG2"]] == 1 | IWS_cluster[["IDH2"]] == 1,  clust_comb1)



## cl9
table(IWS_cluster[["TET2bi"]] == 1,  clust_comb1)

## cl5
table(IWS_cluster[["complex"]] == 1,  clust_comb1)


## GESMD clustering
gesmd_dataset <- gesmd %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("MDS-LB", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(WBC < 20 & ANC < 40)


## Select GESMD variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' - Mutations
kars_gesmd <-  colnames(gesmd_dataset)[43:72]
kars_gesmd <- kars_gesmd[!kars_gesmd %in% c("Otra especificar", "N_aberrations")]
muts_gesmd <- c(colnames(gesmd_dataset)[74:286])

## Remove two mutations not present in any sample
muts_gesmd <- muts_gesmd[!muts_gesmd %in% c("TP53loh", "MLL_PTD", "TP53other")]
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
low_freq_vars_gesmd <- c(low_freq_vars_gesmd, "General Signaling", "Splicing Factors Core")
low_freq_vars_gesmd[is.na(low_freq_vars_gesmd)] <- "del17_17p"
gesmd_dataset_filt <- gesmd_dataset %>%
    select(any_of(c(clin_vars, kars_gesmd, muts_gesmd, gene_cats, "ID"))) %>%
    select(-all_of(c(low_freq_vars_gesmd, "TET2", "TP53", "CUX1", "IDH1", "IDH2", "RET"))) %>%
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
    type = c("s", "n"), graph = FALSE, ncp = 12)

png("figures/GESMD_IWS_clustering/clustering/GESMD_MFA_eigen.png")
plot(GESMD_mfa$global.pca$eig[, 2])
abline(h = 1/length(gesmd_order)*100, col = "red", lty = 2)
dev.off()


GESMD_MFA_sel <- GESMD_mfa$ind$coord[, 1:7]
dist_comb_GESMD_pc1 <- dist(GESMD_MFA_sel, method = "euclidean")
set.seed(27)
sil_scores_gesmd <- sapply(2:20, silhouette_score,
      mat = GESMD_MFA_sel, data_dist = dist_comb_GESMD_pc1)
png("figures/GESMD_IWS_clustering/clustering/GESMD_MFA_kmeans_silhouette.png")
plot(2:20, sil_scores_gesmd, type="b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Silhouette Score",
     main = "Combined clustering")
dev.off()

## Bests clusters: 12
set.seed(27)
clust_gesmd <- kmeans(GESMD_MFA_sel, centers = 8, nstart = 1000)$cluster
tabs_gesmd <- lapply(unlist(avail_vars_list_gesmd[-1]), function(var){
    print(var)
    print(table(GESMD_cluster[[var]], clust_gesmd))
}) 
names(tabs_gesmd) <- unlist(avail_vars_list_gesmd[-1])

## cl6
table(GESMD_cluster[["del5q"]] == 1,  clust_gesmd)
table(GESMD_cluster[["del7"]] == 1,  clust_gesmd)
table(GESMD_cluster[["del5q"]] == 1 & GESMD_cluster[["del7"]] == 1,  clust_gesmd)
table(GESMD_cluster[["del5q"]] == 1 |GESMD_cluster[["del7"]] == 1,  clust_gesmd)

## cl2
table(GESMD_cluster[["Hematopoietic Differentiation"]] == 1,  clust_gesmd)
table(GESMD_cluster[["SF3B1"]] == 1,  clust_gesmd)
table(GESMD_cluster[["SF3B1"]] == 1 |GESMD_cluster[["Hematopoietic Differentiation"]] == 1,  clust_gesmd)

## cl5
table(GESMD_cluster[["TET2bi"]] == 1,  clust_gesmd)
# table(GESMD_cluster[["ZRSR2"]] == 1,  clust_gesmd)
table(GESMD_cluster[["ZRSR2"]] == 1 | GESMD_cluster[["TET2bi"]] == 1,  clust_gesmd)


table(GESMD_cluster[["TP53mono"]] == 1,  clust_gesmd)

# cl4
table(GESMD_cluster[["EZH2"]] == 1,  clust_gesmd)

# cl11
table(GESMD_cluster[["complex"]] == 1,  clust_gesmd)


## cl3
table(GESMD_cluster[["STAG2"]] == 1,  clust_gesmd)

# cl8
table(GESMD_cluster[["BCOR"]] == 1,  clust_gesmd)

# cl10
table(GESMD_cluster[["SRSF2"]] == 1,  clust_gesmd)




IWS_cluster_tree <- IWS_cluster %>%
    mutate(clust = as.factor(clust_comb1))
IWS_tree <- rpart(clust ~ ., data = IWS_cluster_tree, method = "class")

png("figures/GESMD_IWS_clustering/clustering/IWS_tree.png", width = 1000, height = 700)
plot(IWS_tree)
text(IWS_tree, use.n = TRUE)
dev.off()

gesmd_cluster_tree <- GESMD_cluster %>%
    mutate(clust = as.factor(clust_gesmd))
GESMD_tree <- rpart(clust ~ ., data = gesmd_cluster_tree, method = "class")

png("figures/GESMD_IWS_clustering/clustering/GESMD_tree.png", width = 1000, height = 700)
plot(GESMD_tree)
text(GESMD_tree, use.n = TRUE)
dev.off()



classifySamplesIWS <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "complex",
    df$del5q == 1 | df$SF3B1 == 1      ~ "del5q/SF3B1",
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    df$del7 == 1 | df$SETBP1 == 1      ~ "del7/SETBP1",
    df$STAG2 == 1       ~ "STAG2",
    df$`Other Functions` == 1 ~ "Other Functions",
    df$PHF6 == 1        ~ "PHF6",
    TRUE                ~ "Other" 
  )
}
iws_class <- classifySamplesIWS(IWS_cluster)
iws_tab <- table(iws_class, clust_comb1)
matrixStats::rowMaxs(iws_tab) / rowSums(iws_tab)

iws_class_tree <- predict(IWS_tree, IWS_cluster, type = "class")

classifySamplesGESMD <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "complex",
    df$del5q == 1 | df$del7 == 1  ~ "del5q/-7",
    df$EZH2 == 1        ~ "EZH2",
    df$STAG2 == 1       ~ "STAG2",
    df$TET2bi == 1 | df$ZRSR2 == 1     ~ "TET2-bi/ZRSR2",
    df$BCOR == 1 ~ "BCOR",
    df$SRSF2 == 1 ~ "SRSF2",
    df$SF3B1 == 1~ "SF3B1",
    TRUE                ~ "Other" 
  )
}

gesmd_class <- classifySamplesGESMD(GESMD_cluster)
gesmd_tab <- table(gesmd_class, clust_gesmd)
matrixStats::rowMaxs(gesmd_tab) / rowSums(gesmd_tab)
gesmd_class_tree <- predict(GESMD_tree, GESMD_cluster, type = "class")

gesmd_class2 <- classifySamplesIWS(GESMD_cluster)
table(gesmd_class2, clust_gesmd)

iws_class2 <- classifySamplesGESMD(IWS_cluster)
table(iws_class2, clust_comb1)


matrixStats::colMaxs(iws_tab) / colSums(iws_tab)
matrixStats::colMaxs(gesmd_tab) / colSums(gesmd_tab)

table(iws_class, clust_comb1)
#                  clust_comb1
# iws_class           1   2   3   4   5   6   7   8   9  10
#   complex           0   0   3   3  38   7   1   1   0   2
#   del5q/SF3B1      78   5   5   0  18   3   0   4   5  12
#   del7/SETBP1       0   1   1   4   0  52   0  10   1   3
#   EZH2              1   1   5   3   0   4  54   5   6   8
#   Other            39  30   6  44   4   1   8 427  19 158
#   Other Functions   9   1  11   1   0   0   0   2   0   5
#   PHF6              0   0   0  30   0   0   0   0   0   0
#   STAG2            10  75   1   5   0   0  13  10   0   2
#   TET2 bi-allelic   2  11   0   6   0   0   0   2 111   4

table(gesmd_class, clust_gesmd)
#                clust_gesmd
# gesmd_class       1   2   3   4   5   6   7   8   9  10  11  12
#   BCOR            0   4   2   0   0   0   0  14   0   0   0   0
#   complex         0   0   0   0   0   0   0   0   0   0  18   0
#   del5q/-7        0   3   1   1   0  27   3   0   0   0   0   0
#   EZH2            1   1   4  20   2   0   3   0   0   0   0   7
#   Other          26   8   2   2   5   1  21  11  19   5   2 112
#   SF3B1           1  12   0   0   0   0   0   0   0   0   0   0
#   SRSF2           7   0   3   0   0   0   0   1   4  30   0  10
#   STAG2           2   2  28   0   0   0   0   1   0   0   0   2
#   TET2-bi/ZRSR2   4   1   2   2  47   0   0   3   1   3   0   6


classifySamplesComb <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "complex",
    df$del5q == 1      ~ "del5q",
    df$SF3B1 == 1      ~ "SF3B1",
    df$del7 == 1       ~ "-7",
    df$EZH2 == 1        ~ "EZH2",
    df$STAG2 == 1       ~ "STAG2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    TRUE                ~ "Other" 
  )
}

iws_class_comb <- classifySamplesComb(IWS_cluster)
iws_tab_comb <- table(iws_class_comb, clust_comb1)

gesmd_class_comb <- classifySamplesComb(GESMD_cluster)
gesmd_tab_comb <- table(gesmd_class_comb, clust_gesmd)

iws_tab_comb
#                  clust_comb1
# iws_class_comb      1   2   3   4   5   6   7   8   9  10
#   -7                1   0   0   2   0  35   1   8   1   0
#   complex           0   0   3   3  38   7   1   1   0   2
#   del5q            25   0   0   0  16   2   0   4   0   3
#   EZH2              0   1   5   3   0   2  53   5   6   8
#   Other            48  31  18  77   4  18   8 431  19 166
#   SF3B1            53   5   5   0   2   1   0   0   5   9
#   STAG2            10  87   1   5   0   2  13  10   1   2
#   TET2 bi-allelic   2   0   0   6   0   0   0   2 110   4

gesmd_tab_comb
#                  clust_gesmd
# gesmd_class_comb    1   2   3   4   5   6   7   8   9  10  11  12
#   -7                0   1   0   1   0   6   3   0   0   0   0   0
#   complex           0   0   0   0   0   0   0   0   0   0  18   0
#   del5q             0   0   1   0   0  21   0   0   0   0   0   0
#   EZH2              1   0   4  20   2   0   3   0   0   0   0   7
#   Other            35   8   7   4  15   1  21  26  24  36   2 124
#   SF3B1             1  22   0   0   0   0   0   0   0   0   0   0
#   STAG2             2   0  28   0   0   0   0   1   0   0   0   2
#   TET2 bi-allelic   2   0   2   0  37   0   0   3   0   2   0   4

save(sil_scores1, IWS_mfa, clust_comb1, IWS_dataset, IWS_cluster,
     sil_scores_gesmd, GESMD_mfa, clust_gesmd, gesmd_dataset_filt, GESMD_cluster,
    file = "results/clustering/MFA_results.Rdata")

