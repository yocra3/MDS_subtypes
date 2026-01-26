#' ---------------------------
#'
#' Purpose of script:
#'
#'  Test interactions between subgroups of MDS and mutations
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for MDS subgroups interactions
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(cowplot)
library(survminer)
library(survival)
library(tidyverse)

classifySamples <- function(df){
    new_class <- case_when(
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    df$del7 == 1       ~ "7-",
    df$STAG2 == 1       ~ "STAG2",
    df$PHF6 == 1        ~ "PHF6",
    df$BM_BLAST <= 5    ~ "Low blasts",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c( "EZH2", "TET2 bi-allelic", "7-", "STAG2", "PHF6", 
      "Low blasts", "MDS-IB1", "MDS-IB2"))
}

colors6 <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#999999", "grey40",  "black")

## Load data
load("results/preprocess/clinical_preproc.Rdata")
IWS_dataset <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>% ## Filtro IWS
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q 
    filter(if_all(c(BM_BLAST, EZH2, STAG2, del7, TET2bi, PHF6), ~ !is.na(.))) %>%
    mutate(sub_group = classifySamples(.))


IWS_mds <- clinical %>%
    mutate(PLT = log(PLT)) %>% 
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) ## Filtro MDS

## Merge datasets
load("data/GESMD/gesmd_data_1125.Rdata")
gesmd <- gesmd_data_1125 %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

gesmd_dataset <- gesmd %>%
    mutate(PLT = log(PLT)) %>% 
    filter(consensus %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(SF3B1 == 0 & del5q == 0) %>% ## Remove SF3B1 and del5q
    filter(WBC < 20 & ANC < 40) %>%
    filter(if_all(c(BM_BLAST, EZH2, STAG2, del7, TET2bi, PHF6), ~ !is.na(.))) %>%
    mutate(sub_group = classifySamples(.))


gesmd_mds <- gesmd %>%
    mutate(PLT = log(PLT)) %>% 
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>% ## Filtro GESMD
    filter(WBC < 20 & ANC < 40)

## Compute HR for mutations and sub-groups
getEstimates <- function(gene, df){
  mod <- summary(coxph(formula(paste("Surv(OS_YEARS,OS_STATUS) ~", gene, "+ AGE")), 
                       df))
  coef <- mod$coefficients
  conf <- mod$conf.int
  hr <- conf[1, 1]
  hr_l <- conf[1, 3]
  hr_h <- conf[1, 4]
  pval <- coef[1, 5]
  freq <- mean(df[[gene]])
  N <- sum(df[[gene]])
  data.frame(HR = hr, HR_L = hr_l, HR_H = hr_h, Pvalue = pval, 
             Gene = gene, Freq = freq, N = N )
}

## Most frequent genes (>4%)
test_muts <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "CBL", "IDH1")

## Select genes frequent in new sub-groups
test_muts <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1",  "CBL")

groups <- levels(IWS_dataset$sub_group)

mut_hr_clust <- lapply(test_muts, function(gene){
  lapply(groups, function(cl){
    df <- subset(IWS_dataset, sub_group == cl)
    df_est <- getEstimates(gene, df) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_iws <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, IWS_mds) %>%
      mutate(sub_group = "IPSSM cohort")
  }) %>% Reduce(f = rbind)

mut_hr_comb <- rbind(mut_hr_clust, mut_hr_iws) %>%
  mutate(sub_group = factor(sub_group, levels = c(groups, "IPSSM cohort"))) %>%
  filter(N >= 10) %>%
  mutate(sub_group = droplevels(sub_group)) %>%
  as_tibble()

png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations.png", res = 300, heigh = 1300, width = 2200)
mut_hr_comb %>%
  ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-Group") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = c(colors6, "white")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())     
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ STAG2*ASXL1 + AGE, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ plus8*IDH2 + AGE + SEX, IWS_mds))

mut_hr_comb2 <- mut_hr_comb %>%
 group_by(Gene) %>% 
 mutate(HR_L_Ref = HR_L[sub_group == "IPSSM cohort"],
        HR_H_Ref = HR_H[sub_group == "IPSSM cohort"]) %>%
 ungroup()
top_low <- filter(mut_hr_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
    arrange(HR_H - HR_L_Ref)
top_high <- filter(mut_hr_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
    arrange(HR_H_Ref - HR_L)
filter(mut_hr_comb2, HR_H < HR_L_Ref | HR_L > HR_H_Ref)


createSurvFit <- function(df, gene){

    mod <- summary(coxph(formula(paste("Surv(OS_YEARS,OS_STATUS) ~", gene, "+ AGE + SEX")), 
                       df))
    hr <- mod$coefficients[1, 2]
    pval <- mod$coefficients[1, 5]
    survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ gene, df) %>%
        ggsurvplot(data = df, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 pval = sprintf("HR = %.2f\nP = %.3f", hr, pval),
                 pval.coord = c(max(df$OS_YEARS, na.rm = TRUE)*0.7, 0.9),
                 legend.labs  = c("mut", "wt")) +
        xlab("Time (Years)") 
}
makePanelPlot <- function(plot, title){
    plot_grid(plot$plot + 
                theme(legend.position = "none", 
                    plot.title = element_text(hjust = 0.5)) +
                ylab("OS probability") +
                xlab("Time (years)") +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5)), 
                plot$table, ncol = 1, rel_heights = c(1, 0.4))

}

plotPair <- function(group, gene){
    df <- IWS_dataset %>%
        filter(sub_group == group) %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))
    
    p <- createSurvFit(df, gene)
    df2 <- IWS_mds %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))
    p2 <- createSurvFit(df2, gene)

    df3 <- gesmd_dataset %>%
        filter(sub_group == group) %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))   
    p3 <- createSurvFit(df3, gene)

    df4 <- gesmd_mds %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))
    p4 <- createSurvFit(df4, gene)


    plot_grid(
        plot_grid(
            makePanelPlot(p, paste(gene, "in", group, "- IWS")),
            makePanelPlot(p2, paste(gene, "in all IWS")),
            ncol = 2, rel_widths = c(0.8, 1)),
        plot_grid(
            makePanelPlot(p3, paste(gene, "in", group, "- GESMD")),
            makePanelPlot(p4, paste(gene, "in all GESMD")),
            ncol = 2, rel_widths = c(0.8, 1)),
    ncol = 1)
}

# png("figures/GESMD_IWS_clustering/subgroup_interaction/TETbi_ASXL1.png", res = 300, heigh = 3000, width = 2800)
# plotPair("TET2 bi-allelic", "ASXL1")
# dev.off()

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, IWS_dataset, subset = sub_group == "TET2 bi-allelic"))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, IWS_dataset))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, IWS_mds))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1*TET2bi + AGE, IWS_mds))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, gesmd_dataset, subset = sub_group == "TET2 bi-allelic"))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, gesmd_dataset))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE, gesmd_mds))

# png("figures/GESMD_IWS_clustering/subgroup_interaction/TETbi_DNMT3A.png", res = 300, heigh = 1500, width = 2800)
# plotPair("TET2 bi-allelic", "DNMT3A")
# dev.off()

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_dataset, subset = sub_group == "TET2 bi-allelic"))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_dataset))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_mds))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*TET2bi + AGE + SEX, IWS_mds))


# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_dataset, subset = sub_group == "TET2 bi-allelic"))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_dataset))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_mds))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*TET2bi + AGE + SEX, gesmd_mds))

## EZH2 and TET2other
png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_TET2other.png", res = 300, heigh = 3000, width = 2800)
plotPair("EZH2", "TET2other")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other*EZH2 + AGE + SEX, IWS_mds))


summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other*EZH2 + AGE + SEX, gesmd_mds))

## STAG2 and RUNX1
png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_RUNX1.png", res = 300, heigh = 3000, width = 2800)
plotPair("STAG2", "RUNX1")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, IWS_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1*STAG2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, gesmd_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1 + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1*STAG2 + AGE + SEX, gesmd_mds))

## Significativo en IWS pero no en GESMD

## EZH2 and SETBP1
png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_SETBP1.png", res = 300, heigh = 3000, width = 2800)
plotPair("EZH2", "SETBP1")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, IWS_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1*EZH2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, gesmd_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1 + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ SETBP1*EZH2 + AGE + SEX, gesmd_mds))

## STAG2 and U2AF1
png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_U2AF1.png", res = 300, heigh = 3000, width = 2800)
plotPair("STAG2", "U2AF1")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1*STAG2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1*STAG2 + AGE + SEX, gesmd_mds))

## Pocas muestras en GESMD

## EZH2 and U2AF1
png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_U2AF1.png", res = 300, heigh = 3000, width = 2800)
plotPair("EZH2", "U2AF1")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1*EZH2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_dataset, subset = sub_group == "EZH2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ U2AF1*EZH2 + AGE + SEX, gesmd_mds))

## Pocas muestras para ver interacción en GESMD

## STAG2 and ASXL1
png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_ASXL1.png", res = 300, heigh = 3000, width = 2800)
plotPair("STAG2", "ASXL1")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, IWS_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1*STAG2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, gesmd_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1 + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1*STAG2 + AGE + SEX, gesmd_mds))

## La interacción no es significativa

## STAG2 and TET2other
png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_TET2other.png", res = 300, heigh = 3000, width = 2800)
plotPair("STAG2", "TET2other")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other*STAG2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2other*STAG2 + AGE + SEX, gesmd_mds))

## Interacción significativa en GESMD pero marginal en IWS

## STAG2 and DNMT3A
png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_DNMT3A.png", res = 300, heigh = 3000, width = 2800)
plotPair("STAG2", "DNMT3A")
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, IWS_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*STAG2 + AGE + SEX, IWS_mds))

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_dataset, subset = sub_group == "STAG2"))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_dataset))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A + AGE + SEX, gesmd_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*STAG2 + AGE + SEX, gesmd_mds))

## La intersección no es significativa
png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_sel.png", res = 300, heigh = 1300, width = 2200)
mut_hr_comb %>%
    filter(Gene %in% c("TET2other", "RUNX1", "SETBP1", "U2AF1", "ASXL1", "DNMT3A")) %>%
    filter(sub_group %in% c("EZH2", "STAG2", "Low blasts", "IPSSM cohort")) %>%
  ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-Group") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = c(colors6[c(1, 4, 6)], "white")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())     
dev.off()