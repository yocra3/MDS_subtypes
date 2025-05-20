#' ---------------------------
#'
#' Purpose of script:
#'
#'  Exploratory analysis of clustering results using PCA and k-means
#'  clustering. 
#' ---------------------------
#'
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

# Load libraries and data
library(tidyverse)
library(cowplot)
library(pheatmap)
library(survminer)
library(survival)

## Clinical variables
load("results/clustering/lowblast_filter_pc_clustering_11C.Rdata")
load("results/preprocess/clinical_preproc.Rdata")

clinical_low <- clust11_df %>%
    mutate(clust = factor(clust11)) %>%
    left_join(select(clinical, starts_with("OS"), "IPSSM", "ID", "AGE"), by = "ID")

clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")
names(clin_vars) <- clin_vars
low_clin_plots <- lapply(clin_vars, function(var){
    ggplot(clinical_low, aes(x = .data[[var]])) +
        geom_histogram() +
        facet_grid(clust ~ .) +
        theme_bw()
})
png("figures/lowblast_filter_pc_cluster/11_clusters/lowblasts_clinical.png", width = 1000, height = 1000)
plot_grid(plotlist = low_clin_plots, nrow = 3)
dev.off()

png("figures/lowblast_filter_pc_cluster/11_clusters/lowblasts_clinical_violin.png", width = 1000, height = 1000)
clinical_low %>%
    select(all_of(c(clin_vars, "clust"))) %>%
    pivot_longer(cols = !clust, names_to = "Cell", values_to = "value") %>%
    ggplot(aes(y = value, x = clust)) + 
    geom_violin() +
    facet_wrap(~Cell, scales = "free_y") +
    theme_bw()
dev.off()



lapply(clin_vars[clin_vars != "PLT"], function(var){
    summary(glm(formula (paste(var, " ~ clust")), clinical_low, family = "poisson"))
})
summary(lm(PLT ~ clust, clinical_low)) 


## Karyotype events
kar_events <- c("delY", "del11q", "del20q", "del7q",  "plus8",  "del7")
names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(clinical_low[[x]], clinical_low$clust))
kar_tabs_lowblasts <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/lowblast_filter_pc_cluster/11_clusters/lowblasts_karevents_prop.png")
pheatmap(mat = kar_tabs_lowblasts, display_numbers = TRUE)
dev.off()

## Mutations
sel_muts <- c("TET2", "ASXL1", "SRSF2", "DNMT3A", "RUNX1", 
              "STAG2", "U2AF1", "EZH2", "ZRSR2", "TET2bi", "TET2other")
names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(clinical_low[[x]], clinical_low$clust))

mut_tabs_lowblasts <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
png("figures/lowblast_filter_pc_cluster/11_clusters/lowblasts_mutevents_prop.png")
pheatmap(mat = mut_tabs_lowblasts, display_numbers = TRUE)
dev.off()


## IPSSM
png("figures/lowblast_filter_pc_cluster/11_clusters/lowblast_clust_ipssm.png")
table(clinical_low$clust, clinical_low$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/lowblast_filter_pc_cluster/11_clusters/lowblast_clust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_low$clust, clinical_low$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

## Survival analysis
cox_surv_raw <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE , clinical_low )

surv_low <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/lowblast_filter_pc_cluster/11_clusters/lowblasts_survival.png")
surv_low$plot 
dev.off()

cox_aml_raw <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust + AGE + SEX , clinical_low )
cox_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust + AGE + SEX + IPSSM_SCORE, clinical_low )
aml_low <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ clust, clinical_low) %>%
    ggsurvplot(data = clinical_low)

png("figures/lowblast_pc_cluster/lowblasts_amlt.png")
aml_low$plot
dev.off()


## Mutations
clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*clust + AGE + SEX + IPSSM_SCORE, data = .) 
clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*clust + AGE + SEX + IPSSR_SCORE, data = .) 
clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    mutate(clust = relevel(clust, 6)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1*clust + AGE + SEX + IPSSR_SCORE, data = .) 
clinical_low %>%
    mutate(clust = relevel(clust, 6)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(5, 4, 7, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ STAG2*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(5, 4, 7, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ STAG2*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(5, 4, 7, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    mutate(clust = relevel(clust, "4")) %>%
   coxph(Surv(OS_YEARS,OS_STATUS) ~ EZH2*clust + AGE + SEX + IPSSR_SCORE, data = .) 
clinical_low %>%
    filter(clust %in% c(5, 4, 7, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    mutate(clust = relevel(clust, "4")) %>%
   coxph(Surv(OS_YEARS,OS_STATUS) ~ EZH2*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    filter(clust != 6) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*clust + AGE + SEX + IPSSR_SCORE, data = .) 
clinical_low %>%
    filter(clust != 6) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(5, 6, 2, 7)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ZRSR2*clust + AGE + SEX + IPSSR_SCORE, data = .) 

clinical_low %>%
    filter(clust %in% c(5, 6, 2, 7)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ZRSR2*clust + AGE + SEX + IPSSM_SCORE, data = .) 

clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*clust + AGE + SEX + IPSSR_SCORE, data = .) 
clinical_low %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*clust + AGE + SEX + IPSSM_SCORE, data = .) 
