#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define clusters of patients using FactoMineR. Not include individual karyotypes
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

## Select variables for clustering
#' - Demographic: age 
#' - Clinical: Blood cell proportions
#' - Karyotype events: selected events and number of aberrations
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clust_dataset <- clinical_all %>%
    select(c("AGE", "SEX", clin_vars))



## Run Multiple Factor Analysis (MFA) with FactoMiner
mfa_nokar <- MFA(clust_dataset, group = c(1, 1, length(clin_vars)),
    type = c("s", "n", "s"), ncp = 8)

## Make clusters with Hierarchical Clustering on Principal Components (HCPC)
mfa_clust_nokar <- HCPC(mfa_nokar,  nb.clust = 6)

png("inertia_nokar.png")
plot(mfa_clust_nokar$call$t$inert.gain[1:50])
dev.off()

load("results/clustering/factominer_categorical.Rdata")

clust <- mfa_clust$data.clust$clust

clust_nokar <- mfa_clust_nokar$data.clust$clust

table(clust, clust_nokar)




## Compare clustering with other variables
table(clust_nokar, clinical_all$SEX)
chisq.test(table(clust, clinical_all$SEX))

table(clust_nokar, clinical_all$IPSSM)
chisq.test(table(clust, clinical_all$IPSSM))
kar_events <- c("monoy", "del11q", "del5q", "del12p",
                "del20q", "del7q", "triso8", "triso19", "mono7", 
                "i17q", "chr3ab")

names(kar_events) <- kar_events
lapply(kar_events, function(x) table(clinical_all[[x]], clust))

### Combine clusters 8 and 10 - "controls"
clin_clust <- mutate(clinical_all, 
    kar_clust = factor(clust, levels = c(1:16)), 
    clust_nokar = factor(clust_nokar, levels = c(8, 1:7, 9:10)))
coxph(Surv(OS_YEARS,OS_STATUS) ~ kar_clust*hma + IPSSM_SCORE, clin_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust_nokar*hma + IPSSM_SCORE, clin_clust)


save(mfa, mfa_clust, file = "results/clustering/factominer_categorical.Rdata")