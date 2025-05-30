#' ---------------------------
#'
#' Purpose of script:
#'
#'  Descriptives of the cohort
#' 
#' ---------------------------
#'
#' Notes:
#' Compute the descriptives of GESMD and IWS cohorts. Include patients with
#' classification and patients used for clustering. Also classify patients
#' in the larger cohorts.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(tidyverse)

load("results/clustering/combined_pc_clustering.Rdata")

## Compute descriptives
getIQR <- function(vec){
   sprintf("%.1f (%.1f-%.1f)", 
           median(vec, na.rm = TRUE), 
           quantile(vec, probs = 0.25, na.rm = TRUE),
           quantile(vec, probs = 0.75, na.rm = TRUE))
}


getProp <- function(cond, var){
    sprintf("%i (%.1f%%)", sum(cond, na.rm = TRUE), mean(!is.na(var) & cond)*100)
}

getPropNA <- function(var){
    sprintf("%i (%.1f%%)", sum(is.na(var)), mean(is.na(var))*100)
}

summarize_fun <- function(df){
     summarize(df, 
            N = n(),
            Females = getProp(SEX == "F", SEX),
            Males = getProp(SEX == "M", SEX),
            Age = getIQR(AGE),
            `BM Blasts` = getIQR(BM_BLAST),
            `WBC count` = getIQR(WBC),
            `Neutrophil Count` = getIQR(ANC),
            `Monocyte Count` = getIQR(MONOCYTES),
            HB = getIQR(HB),
            PLT = getIQR(PLT),
            `Low blasts` = getProp(consensus == "Low blasts", consensus),
            `MDS-IB1` = getProp(consensus == "MDS-IB1", consensus),
            `MDS-IB2` = getProp(consensus == "MDS-IB2", consensus),
            `IPSSM Very-Low` = getProp(IPSSM == "Very-Low", IPSSM),
            `IPSSM Low` = getProp(IPSSM == "Low", IPSSM),
            `IPSSM Moderate-Low` = getProp(IPSSM == "Moderate-Low", IPSSM),
            `IPSSM Moderate-High` = getProp(IPSSM == "Moderate-High", IPSSM),
            `IPSSM High` = getProp(IPSSM == "High", IPSSM),
            `IPSSM Very-High` = getProp(IPSSM == "Very-High",  IPSSM),
            IPSSM_NA = getPropNA(IPSSM)
           )
}


descriptives_cluster <- combined_full %>%
    mutate(IPSSM = ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_ALTO", "High",
        ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_BAJO", "Low", 
           ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_ALTO", "Moderate-High",
            ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_BAJO", "Moderate-Low", IPSSM)) )))  %>%
    group_by(dataset) %>%
    summarize_fun() %>% 
    t()

write.table(descriptives_cluster, 
            file = "results/GESMD_IWS_clustering/clust_samples_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)

## Tests
poisson_test <- function(var, df){
    poiss_lm <- summary(glm(formula (paste(var, " ~ dataset")), df, 
                      family = "poisson"))
    poiss_lm$coefficients[2, 4]
}

lm_test <- function(var, df){
    lm_res <- summary(lm(formula (paste(var, " ~ dataset")), df))
    lm_res$coefficients[2, 4]
}

chisq_test <- function(var, df){
    chisq_res <- chisq.test(table(df[[var]], df$dataset))
    chisq_res$p.value
}
# Perform tests
clust_test <- c(sapply(c("SEX", "consensus", "IPSSM"), chisq_test, df = combined_full),
            sapply(c("AGE", "HB", "PLT"), lm_test, df = combined_full),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = combined_full))
clust_test


## Compute descriptives for GESMD and IWS full cohorts
load("results/preprocess/clinical_preproc.Rdata")
load("data/GESMD/gesmd_data_all.Rdata")

classify_patient <- function(df){
    classification <- ifelse(df$delY == 1, "Y-",
        ifelse(df$del7 == 1, "7-",
           ifelse(df$complex == 1, "complex",
              ifelse(df$del20q == 1, "del20q",
                ifelse(df$plus8 == 1, "8+", 
                    ifelse(df$del7q == 1, "del7q",
                        ifelse(df$TET2bi == 1, "TET2 bi-allelic",
                            ifelse(df$TET2other == 1, "TET2 monoallelic",
                                ifelse(df$STAG2 == 1, "STAG2",
                                    ifelse(df$WBC > 6, "Midly Leukopenic",                 
                                           "Highly Leukopenic"))))))))))
    classification
}

## Merge datasets
gesmd <- gesmd_data %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017) %>%
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>%
    filter(SF3B1 == 0 & del5q == 0) %>%
    filter(WBC < 20) %>%
    mutate(sub_group = classify_patient(.)) %>%
    filter(!is.na(sub_group)) %>%
    mutate(IPSSM = ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_ALTO", "High",
        ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_BAJO", "Low", 
           ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_ALTO", "Moderate-High",
            ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_BAJO", "Moderate-Low", 
                ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MUY_BAJO", "Very-Low", 
                    ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MUY_ALTO", "Very-High", IPSSM)) )))))



clinical_blasts <- clinical %>%
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>%
    filter(SF3B1 == 0 & del5q == 0) %>%
    mutate(complex = ifelse(complex == "complex", 1, 0)) %>%
    mutate(sub_group = classify_patient(.)) %>%
    filter(!is.na(sub_group)) 

## Combine datasets for descriptives
full_dataset <- rbind(gesmd %>% mutate(dataset = "GESMD") %>%
    select(names(clust_test), dataset), 
    clinical_blasts %>% mutate(dataset = "IWS") %>%
    select(names(clust_test), dataset))

descriptives_full <- full_dataset %>%
    group_by(dataset) %>%
    summarize_fun() %>% 
    t()

write.table(descriptives_full, 
            file = "results/GESMD_IWS_clustering/all_samples_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)

all_test <- c(sapply(c("SEX", "consensus", "IPSSM"), chisq_test, df = full_dataset),
            sapply(c("AGE", "HB", "PLT"), lm_test, df = full_dataset),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = full_dataset))
all_test

save(gesmd, clinical_blasts, file = "results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")