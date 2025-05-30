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
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.0 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(FactoMineR)
library(survival)
library(survminer)

load("results/preprocess/clinical_preproc.Rdata")
load("results/clustering/clinical_clustering.Rdata")
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
    mutate(Sex = factor(SEX)) %>%
    select(c(Sex, AGE, clin_vars, kar_events, "N_aberrations"))
clust_dataset[kar_events] <- lapply(kar_events, function(x){
   vec = clust_dataset[[x]]
   factor(c("no", "yes")[vec + 1])
})

clust_dataset <- clinical_all %>%
    mutate(Sex = factor(SEX)) %>%
    select(c(Sex, AGE, clin_vars, "N_aberrations"))

## Run Multiple Factor Analysis (MFA) with FactoMiner
mfa <- MFA(clust_dataset, group = c(1, 1, lengths(list(clin_vars, kar_events)), 1),
    type = c("n", "s", "s", "n", "s"))

mfa <- MFA(clust_dataset, group = c(2, lengths(list(clin_vars)), 1),
    type = c("m", "s", "s"))


png("figures/all_clust_eigen.png")
plot(mfa$global.pca$eig[, 3])
dev.off()


clust_dataset2 <- clinical_all %>%
    mutate(Sex = ifelse(SEX == "F", 0, 1)) %>%
    select(c(Sex, AGE, clin_vars, kar_events, "N_aberrations"))
pca_all <- PCA(clust_dataset2)
pca_clust <- HCPC(pca_all,  nb.clust = 11)

## Make clusters with Hierarchical Clustering on Principal Components (HCPC)
mfa_clust <- HCPC(mfa,  nb.clust = 9)

png("figures/all_clust_inertia.png")
plot(mfa_clust$call$t$inert.gain[1:30])
dev.off()

png("figures/all_pca_clust_inertia.png")
plot(pca_clust$call$t$inert.gain[1:30])
dev.off()

clust <- mfa_clust$data.clust$clust
clust <- pca_clust$data.clust$clust

## Compare clustering with other variables
table(clust, clin_clust)
chisq.test(table(clust, clin_clust))


table(clust, clinical_all$SEX)
chisq.test(table(clust, clinical_all$SEX))

table(clust, clinical_all$IPSSM)
chisq.test(table(clust, clinical_all$IPSSM))

names(kar_events) <- kar_events
lapply(kar_events, function(x) table(clinical_all[[x]], clin_clust))
lapply(kar_events, function(x) table(clinical_all[[x]], clust))


## Comparison with survival
clinical_clust <- mutate(clinical_all, all_clust = clust, clin_clust = clin_clust)
clinical_clust$clust <- relevel(clinical_clust$all_clust, "8")
png("Surv_clust.png")
os_clust <- survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()

coxph(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*hma + IPSSM_SCORE, clinical_clust)
## Different clusters: 1, 3, 6
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE*hma, clinical_clust)

png("figures/Surv_clust6_hma.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == "6" ), data = clinical_clust) 
dev.off()

png("figures/Surv_hma.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust), data = clinical_clust) 
dev.off()

png("figures/Surv_clust1_3_6_hma.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust %in% c("1", "3", "6")), data = clinical_clust) 
dev.off()

clust_groups <- c("1", "3", "6")
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust , subset = clust %in% clust_groups)
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + IPSSM_SCORE + AGE + SEX, clinical_clust , subset = clust %in% clust_groups)

coxph(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + IPSSM_SCORE + AGE + SEX, clinical_clust)

clin_clust2 <- mutate(clinical_clust, hma_clust = ifelse(clust %in% clust_groups, "hma_clust", "Other"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*hma_clust + IPSSM_SCORE + AGE + SEX, clin_clust2)

coxph(Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = IPSSM %in% c("High", "Very-High"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma, clin_clust2, subset = IPSSM %in% c("High", "Very-High") & hma_clust == "hma_clust")
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma, clin_clust2, subset = IPSSM %in% c("High", "Very-High") & hma_clust == "Other")

coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + IPSSM + AGE + SEX , clin_clust2, subset = IPSSM %in% c("High", "Very-High") & hma_clust == "hma_clust")
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + IPSSM + AGE + SEX, clin_clust2, subset = IPSSM %in% c("High", "Very-High") & hma_clust == "Other")

clin_clust2$hmaclust_int <- paste(clin_clust2$hma_clust, clin_clust2$hma )
clin_clust2$hmaclust_int2 <- paste(clin_clust2$clust == "6", clin_clust2$hma )

png("figures/Surv_hmaclust_hma.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hmaclust_int, clin_clust2), data = clin_clust2) 
dev.off()


png("figures/Surv_hmaclust_hma_highrisk.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hmaclust_int, clin_clust2,  subset = IPSSM %in% c("High", "Very-High")), data = clin_clust2) 
dev.off()

png("figures/Surv_clust6_hma_highrisk.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hmaclust_int2, clin_clust2,  subset = IPSSM %in% c("High", "Very-High") & (hma_clust == "Other" | clust == "6")), data = clin_clust2) 
dev.off()

clin_clust2$c3_int <- paste(clin_clust2$clust == "3", clin_clust2$hma )

png("figures/Surv_clust3_hma_highrisk.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ c3_int, clin_clust2,  subset = IPSSM %in% c("High", "Very-High") & (hma_clust == "Other" | clust == "3")), data = clin_clust2) 
dev.off()

clin_clust2$c1_int <- paste(clin_clust2$clust == "1", clin_clust2$hma )

png("figures/Surv_clust1_hma_highrisk.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ c1_int, clin_clust2,  subset = IPSSM %in% c("High", "Very-High") & (hma_clust == "Other" | clust == "1")), data = clin_clust2) 
dev.off()




png("figures/Surv_hmaclust_hmatreated_high.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hma_clust, clin_clust2,  subset = IPSSM %in% c("High", "Very-High") & hma == 1), data = clin_clust2) 
dev.off()


png("figures/Surv_clust_hmatreated_high.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clin_clust2,  subset = IPSSM %in% c("High", "Very-High") & hma == 1), data = clin_clust2) 
dev.off()

coxph(Surv(OS_YEARS,OS_STATUS) ~ clust, clin_clust2, subset = IPSSM %in% c("High", "Very-High") & hma == 1)

png("figures/Surv_hmaclust_hma_veryhighrisk.png")
ggsurvplot( survfit(Surv(OS_YEARS,OS_STATUS) ~ hmaclust_int, clin_clust2,  subset = IPSSM == "Very-High"), data = clin_clust2) 
dev.off()

save(pca_clust, pca_all, file = "results/clustering/factominer_pca_all_continuous.Rdata")