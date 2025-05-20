#' ---------------------------
#'
#' Purpose of script:
#'
#'  Subdivide MDS blasts based on interaction effects on survival
#' 
#' ---------------------------
#'
#' Notes:

#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(MASS)
library(cluster)
library(tidyverse)
library(rpart)
library(survival)
library(corrplot)
library(pheatmap)


load("results/preprocess/clinical_preproc.Rdata")

mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")

## Select variables for clustering
#' - Clinical: Blood cell proportions
#' - Karyotype events
#' . Mutations 
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
kar_events <- c("delY", "del11q", "del5q", "del12p", 
                "del20q", "del7q", "plus8", "plus19", "del7")
mutations <- colnames(clinical)[colnames(clinical) %in% colnames(mutation)]
mut_vars <- mutations[mutations != "ID"]

# Select data for computing interaction effects
int_dataset <- clinical %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    select(all_of(c("ID", clin_vars, kar_events, mut_vars, "AGE", "SEX","BM_BLAST")), starts_with("OS")) %>%
    filter(complete.cases(.))  %>%
    mutate(AGE = (AGE - min(AGE))/diff(range(AGE)),
    BM_BLAST = (BM_BLAST - min(BM_BLAST))/diff(range(BM_BLAST))
    )

freq_vars <- select(clinical, all_of(c(kar_events, mut_vars))) %>%
    as.matrix() %>%
    colMeans()

sel_vars <- names(freq_vars)[freq_vars > 0.01]
sel_vars <- sel_vars[!sel_vars %in% c("TET2", "TP53")]

int_coef_df <- lapply(sel_vars, function(var1){
    res <- lapply(c(sel_vars, "AGE", "SEX", "BM_BLAST"), function(var2){
        if (var2 %in% sel_vars & (sum(int_dataset[[var1]] == 1 &  int_dataset[[var2]] == 1) < 10)){
            int_coef <- NA
            var1_coef <- NA
            var2_coef <- NA
            int_p <- NA
        } else {
        surv_mod <- summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ get(var1) * get(var2), 
                data = int_dataset))
        var1_coef <- surv_mod$coefficients[1,1]
        var2_coef <- surv_mod$coefficients[2,1]
          int_coef <- surv_mod$coefficients[3,1]
          int_p <- surv_mod$coefficients[3,5]
        }
        tibble(var1 = var1, var2 = var2, var1_coef = var1_coef, var2_coef = var2_coef, 
            int_coef = int_coef, int_p = int_p) 
    })
    Reduce(rbind, res)
}) %>% Reduce(rbind, .)

int_coef_mat <- int_coef_df %>%
    select(var1, var2, int_coef) %>%
    pivot_wider(names_from = var2, values_from = int_coef) %>%
    column_to_rownames("var1")
int_coef_mat <- int_coef_mat[rowSums(!is.na(int_coef_mat)) > 5, colSums(!is.na(int_coef_mat)) > 6] %>%
    as.matrix()
png("figures/interaction_coef_matrix.png", width = 600, height = 600)
pheatmap(int_coef_mat, display_numbers = TRUE, 
         main = "Interaction coefficients", cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

var1_coef_mat <- int_coef_df %>%
    select(var1, var2, var1_coef) %>%
    pivot_wider(names_from = var2, values_from = var1_coef) %>%
    column_to_rownames("var1")
var1_coef_mat <- var1_coef_mat[rownames(int_coef_mat), colnames(int_coef_mat)]
png("figures/var1_coef_matrix.png", width = 600, height = 600)
pheatmap(var1_coef_mat, display_numbers = TRUE, 
         main = "Row coefficients", cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()


# var2_coef_mat <- int_coef_df %>%
#     select(var1, var2, var2_coef) %>%
#     pivot_wider(names_from = var2, values_from = var2_coef) %>%
#     column_to_rownames("var1")
# var2_coef_mat <- var2_coef_mat[rownames(int_coef_mat), colnames(int_coef_mat)]
# png("figures/var2_coef_matrix.png", width = 600, height = 600)
# pheatmap(var2_coef_mat, display_numbers = TRUE, 
#          main = "Column coefficients", cluster_rows = TRUE, cluster_cols = TRUE)
# dev.off()

top_coefs <- filter(int_coef_df, int_p < 0.05)
sort(table(top_coefs$var1))
    #   CBL   CSNK1A1    del20q     del7q      delY    DNMT3A      EZH2      GNAS 
    #     1         1         1         1         1         1         1         1 
    #  GNB1      IDH2     plus8     PPM1D    PTPN11    SETBP1    TET2bi       WT1 
    #     1         1         1         1         1         1         1         1 
    #  BCOR      NPM1      NRAS     SRSF2     STAG2   MLL_PTD     RUNX1 TET2other 
    #     2         2         2         2         2         3         3         3 
    # U2AF1     ASXL1 
    #     3         5 
sort(matrixStats::rowMads(as.matrix(var1_coef_mat), na.rm = TRUE))
#     TET2bi      RUNX1      ZRSR2      plus8  TET2other      ASXL1     SETBP1 
# 0.03694791 0.04101541 0.04379471 0.04794190 0.07568596 0.07828757 0.07980609 
#     DNMT3A       EZH2      U2AF1      STAG2       PHF6     del20q      SRSF2 
# 0.08342120 0.08438815 0.08568765 0.09067102 0.10192490 0.10379163 0.11229620 
#       del7       BCOR       CUX1      CEBPA       IDH2       IDH1       NRAS 
# 0.13194160 0.15602104 0.16665143 0.21873556 0.22230416 0.22236396 0.24422943 
#       ETV6    MLL_PTD      ETNK1 
# 0.32925338 0.42328113 0.44899262 


## Most interesting genes: SRSF2, STAG2, RUNX1, U2AF1, ASXL1, TET2other
comb_mat <- cbind(int_coef_mat, var1_coef_mat)

penalized_dist <- function(data, penalty = 1.0) {
  n <- nrow(data)
  dmat <- matrix(NA, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      xi <- data[i, ]
      xj <- data[j, ]
      
      both_obs <- !is.na(xi) & !is.na(xj)
      n_obs <- sum(both_obs)
      n_miss <- length(xi) - n_obs
      
      if (n_obs > 0) {
        d <- sqrt(sum((xi[both_obs] - xj[both_obs])^2) + penalty * n_miss)
      } else {
        d <- penalty * length(xi)  # all missing, maximum penalty
      }
      dmat[i, j] <- dmat[j, i] <- d
    }
  }
  as.dist(dmat)
}

penalized_dist2 <- function(data, penalty = 1.0) {
  n <- nrow(data)
  dlist <- list()

  
  for (i in 1:(n-1)) {
    dlist[[i]] <- list()
    for (j in (i+1):n) {
      xi <- data[i, ]
      xj <- data[j, ]
      
      both_obs <- !is.na(xi) & !is.na(xj)
        n_obs <- sum(both_obs)

      if (n_obs > 0) {
        d <- sqrt((xi[both_obs] - xj[both_obs])^2)
      } else {
        d <- NA
      }
      dlist[[i]][[j]] <- d
    }
  }
  dlist
}

# Apply to your data
hclust_coef <- hclust(penalized_dist(int_coef_mat, penalty = 1.5), method = "ward.D2")
#hclust_coef <- hclust(dist(int_coef_mat), method = "ward.D2")

png("figures/coef_hclust.png", width = 600, height = 600)
plot(hclust_coef, labels = rownames(comb_mat))
dev.off()