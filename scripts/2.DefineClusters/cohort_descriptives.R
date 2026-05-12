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
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.8 R
#'
#' ---------------------------


# Load libraries and data
library(tidyverse)
library(ipssm)
library(readxl)

load("results/clustering/MFA_results.Rdata")

classifySamples <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "Complex",
    df$del5q == 1      ~ "del5q-IB",
    df$SF3B1 == 1      ~ "SF3B1-IB",
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2-bi",
    df$del7 == 1       ~ "-7",
    df$STAG2 == 1       ~ "STAG2",
    df$BM_BLAST <= 5    ~ "MDS-LB",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c("EZH2", "TET2-bi",  "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex",
      "MDS-LB", "MDS-IB1", "MDS-IB2"))
}

classifySamplesTaxonomy <- function(df){
    new_class <- case_when(
    df$DDX41 == 1       ~ "DDX41",
    df$NPM1 == 1 | (df$FLT3 == 1 & df$NPM1 == 1) ~ "AML-like",
    df$TP53multi == 1 | df$complex == 1 ~ "TP53-complex",
    df$del7 == 1 | df$SETBP1 == 1       ~ "7-/SETBP1",
    df$del5q == 1      ~ "del(5q)",
    df$EZH2 == 1 & df$ASXL1 == 1  ~ "EZH2-ASXL1",
    df$IDH2 == 1 | df$IDH1 == 1 | (df$STAG2 == 1 & df$ASXL1 == 1) | (df$STAG2 == 1 & df$SRSF2 == 1) ~ "IDH-STAG2",
    df$BCOR == 1 | df$BCORL1 == 1 ~ "BCOR/L1",
    df$TET2bi == 1 | (df$TET2other == 1 & df$SRSF2 == 1) ~ "bi-TET2",
    df$U2AF1 == 1 ~ "U2AF1",
    df$SRSF2 == 1 ~ "SRSF2",
    df$ZRSR2 == 1 ~ "ZRSR2",
    df$SF3B1 == 1 ~ "SF3B1",
    df$DNMT3A == 1 | df$TET2 == 1 | df$TP53mono == 1 | df$delY == 1 ~ "CCUS-like",
    TRUE                ~ "Other" 
  )
  factor(new_class)
}


molecular_class <- read_xlsx("data/BLOOD_BLD-2023-023727-mmc1.xlsx", sheet = "Table S3 (full database)") 




## Preprocess IWS
load("results/preprocess/clinical_preproc.Rdata")
IWS_full <- clinical %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>%
    mutate(PLT2 = pmin(PLT, 250))

IWS_mds <- IWS_full %>%
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(if_all(c(BM_BLAST, EZH2, STAG2, del7, TET2bi, del5q, SF3B1), ~ !is.na(.))) %>%
    mutate(complex = ifelse(complex == "complex", 1, 0)) %>%
    mutate(sub_group = classifySamples(.),
    mol_manual = classifySamplesTaxonomy(.)) %>%
    left_join(molecular_class %>% select(ID, MOLECULAR_GROUP), by = "ID")
IWS_dataset_filt <- subset(IWS_mds, ID %in% IWS_dataset$ID)


## Preprocess GESMD 
load("data/GESMD/gesmd_data_1125.Rdata")
gesmd <- gesmd_data_1125_filt %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

gesmd_full <- gesmd %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>%
    filter(!is.na(TP53)) %>%
    filter(WBC < 20 & ANC < 40) %>%
    mutate(TP53mono = ifelse(TP53 == 1 & TP53multi == 0, 1, 0))  

gesmd_dataset <- gesmd_full %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(if_all(c(BM_BLAST, EZH2, STAG2, del7, TET2bi, del5q, SF3B1, complex), ~ !is.na(.))) %>%
    mutate(sub_group = classifySamples(.),
    mol_manual = classifySamplesTaxonomy(.)) 

ipssm_process <- IPSSMprocess(gesmd_dataset)
ipssm_res <- IPSSMmain(ipssm_process)
ipssm_annot <- IPSSMannotate(ipssm_res)

gesmd_ipssm <-  ipssm_annot %>% 
    mutate(IPSSM = IPSSMcat_mean,
    IPSSM = gsub(" ", "-", IPSSM , fixed = TRUE),
    IPSSM_SCORE = IPSSMscore_mean) %>%
    select(ID, IPSSM, IPSSM_SCORE)

gesmd_dataset <- gesmd_dataset %>%
    rows_patch(y = select(gesmd_ipssm, -IPSSM_SCORE), by = "ID") %>%
    left_join(select(gesmd_ipssm, -IPSSM), by = "ID")

gesmd_cluster <- subset(gesmd_dataset, ID %in% gesmd_dataset_filt$ID)

## MLL
load("results/hershberger/hershberger_mds.Rdata")
load("results/hershberger/hershberger_full.Rdata")

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
            `MDS-LB` = getProp(consensus == "Low blasts", consensus),
            `MDS-IB1` = getProp(consensus == "MDS-IB1", consensus),
            `MDS-IB2` = getProp(consensus == "MDS-IB2", consensus),
            `del5q` = getProp(consensus == "del5q", consensus),
            `SF3B1` = getProp(consensus == "mutated SF3B1", consensus),
            `TP53` = getProp(consensus == "mutated TP53", consensus),
            `IPSSM Very-Low` = getProp(IPSSM == "Very-Low", IPSSM),
            `IPSSM Low` = getProp(IPSSM == "Low", IPSSM),
            `IPSSM Moderate-Low` = getProp(IPSSM == "Moderate-Low", IPSSM),
            `IPSSM Moderate-High` = getProp(IPSSM == "Moderate-High", IPSSM),
            `IPSSM High` = getProp(IPSSM == "High", IPSSM),
            `IPSSM Very-High` = getProp(IPSSM == "Very-High",  IPSSM),
            IPSSM_NA = getPropNA(IPSSM)
           )
}


summarize_fun_input <- function(df, clin_vars, mutations){
        summarize(df, 
        across(all_of(clin_vars),
            ~ getIQR(.),
            .names = "{col}"),
        across(all_of(mutations),
            ~ getProp(. == 1, .),
            .names = "{col}"),
    )  
}

## Union of IWS and GESMD variables
cluster_comb <- rbind(
    IWS_dataset_filt %>% mutate(dataset = "IWS") %>%
    select(SEX, AGE, BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT, consensus, IPSSM, dataset), 
    gesmd_cluster %>% mutate(dataset = "GESMD") %>%
    select(SEX, AGE, BM_BLAST, WBC, ANC, MONOCYTES, HB, PLT, consensus, IPSSM, dataset)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")))

descriptives_cluster <- cluster_comb %>%
     group_by(dataset) %>%
    summarize_fun() %>% 
    t()

write.table(descriptives_cluster, 
            file = "results/GESMD_IWS_clustering/clust_samples_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)


all_variables <- union(colnames(GESMD_cluster), colnames(IWS_cluster))
cluster_comb2 <- bind_rows(
    IWS_cluster %>% mutate(dataset = "IWS"), 
    GESMD_cluster %>% mutate(dataset = "GESMD") 
) %>%
    select(all_of(all_variables), dataset) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")))

clin_vars <- c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
descriptives_cluster2 <- cluster_comb2 %>%
     group_by(dataset) %>%
    summarize_fun_input(clin_vars, all_variables[!all_variables %in% clin_vars]) %>% 
    t() %>%
    as_tibble(rownames = "Variable")


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
clust_test <- c(sapply(c("SEX", "consensus", "IPSSM"), chisq_test, df = cluster_comb),
            sapply(c("AGE", "HB", "PLT"), lm_test, df = cluster_comb),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = cluster_comb))
clust_test

mutations_test <- sapply(all_variables[!all_variables %in% clin_vars], chisq_test, df = cluster_comb2)
descriptives_cluster2$Test <- c(NA, clust_test[clin_vars], mutations_test)

write.table(descriptives_cluster2, 
            file = "results/GESMD_IWS_clustering/clust_samples_input_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)



## Compute descriptives for GESMD and IWS full cohorts
mds_morph_dataset <- bind_rows(
    gesmd_dataset %>% mutate(dataset = "GESMD"),
    IWS_mds %>% mutate(dataset = "IWS"),
    hersh_mds %>% mutate(dataset = "MLL", 
        consensus = ifelse(WHO == "MDS-LB", "Low blasts", WHO),
        IPSSM = gsub(" ", "-", IPSSM))) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD", "MLL")))

descriptives_mds <- mds_morph_dataset %>%
    group_by(dataset) %>%
    summarize_fun() %>% 
    t()

write.table(descriptives_mds, 
            file = "results/GESMD_IWS_clustering/mds_morph_samples_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)

poisson_test_f <- function(var, df){
    poiss_lm <- glm(formula (paste(var, " ~ dataset")), df, 
                    family = "poisson")
    anova(a, test = "Chisq")$`Pr(>Chi)`[2]
}

lm_test_f <- function(var, df){
    lm_res <- summary(lm(formula (paste(var, " ~ dataset")), df))
    pf(lm_res$fstatistic[1], lm_res$fstatistic[2], lm_res$fstatistic[3], lower.tail = FALSE)
}
mds_test <- c(sapply(c("SEX", "consensus", "IPSSM"), chisq_test, df = mds_morph_dataset),
            sapply(c("AGE", "HB", "PLT"), lm_test_f, df = mds_morph_dataset),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = mds_morph_dataset))
mds_test

save(gesmd_dataset, IWS_mds, file = "results/GESMD_IWS_clustering/gesmd_IWS_mds.Rdata")
save(gesmd_full, IWS_full, file = "results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

