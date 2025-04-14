#' ---------------------------
#'
#' Purpose of script:
#'
#' Correlate genetic scores with overall risk
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.4 R
#'
#' ---------------------------
library(survival)
library(tidyverse)
library(corrplot)

load("results/tsne_perp10.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(clust10_tree, levels = 1:17))


bol_scores <- read_table("results/mutations/boolean_scores.tsv")
vaf_scores <- read_table("results/mutations/vaf_scores.tsv")

bol_cor <-  cor(data.matrix(bol_scores[, -60]))
png("figures/NN_models/boolean_scores_cor.png")
corrplot(bol_cor, order = "AOE")
dev.off()


vaf_cor <-  cor(data.matrix(vaf_scores[, -60]))
png("figures/NN_models/vaf_scores_cor.png")
corrplot(vaf_cor, order = "AOE")
dev.off()

## Very high correlation between scores
### Select two scores representatives of the clusters
sel_gos <- c("GO:0051321", "GO:1900225")
load("results/preprocess/clinical_preproc.Rdata")

bol <- left_join(bol_scores, select(clinical_clust, ID, OS_YEARS, OS_STATUS, AGE, SEX, IPSSM_SCORE), by = "ID") %>%
    select(-ID)

coxph(Surv(OS_YEARS, OS_STATUS) ~ `GO:0051321` + `GO:1900225` + SEX + AGE, bol)
coxph(Surv(OS_YEARS, OS_STATUS) ~ `GO:0051321` + `GO:1900225` + SEX + AGE + IPSSM_SCORE, bol)

cor(bol$IPSSM_SCORE, bol$`GO:0051321`, use = "complete")
cor(bol$IPSSM_SCORE, bol$`GO:1900225`, use = "complete")

png("figures/NN_models/boolean_scores_IPSSM.png")
par(mfrow = c(1, 2))
plot(bol$IPSSM_SCORE, bol$`GO:0051321`)
plot(bol$IPSSM_SCORE, bol$`GO:1900225`)
dev.off()

vaf <- left_join(vaf_scores, select(clinical_clust, ID, OS_YEARS, OS_STATUS, AGE, SEX, IPSSM_SCORE), by = "ID") %>%
    select(-ID)

coxph(Surv(OS_YEARS, OS_STATUS) ~ `GO:0051321` + SEX + AGE, vaf)
coxph(Surv(OS_YEARS, OS_STATUS) ~ `GO:0051321` + SEX + AGE + IPSSM_SCORE, vaf)
cor(vaf$IPSSM_SCORE, vaf$`GO:0051321`, use = "complete")

vaf_scores2 <- read_table("results/mutations/vaf_scores2.tsv")

vaf2_cor <-  cor(data.matrix(vaf_scores2[, -11]))
png("figures/NN_models/vaf_scores2_cor.png")
corrplot(vaf2_cor, order = "AOE")
dev.off()

bol_nn <- read_table("results/mutations/boolean_scores_full.tsv")
bol_nn_cor <-  cor(data.matrix(bol_nn[, -60]))
png("figures/NN_models/boolean_nn_cor.png")
corrplot(bol_nn_cor, order = "AOE")
dev.off()
