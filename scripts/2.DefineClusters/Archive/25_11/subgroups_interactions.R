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

load( "results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/preprocess/clinical_preproc.Rdata")

## Prepare data

clinical_blasts <- clinical_blasts %>%
    mutate(sub_group = factor(sub_group,
        levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", 
            "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2")))

clinical_mds <- clinical %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) 

groups <- levels(clinical_blasts$sub_group)
names(groups) <- groups

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

test_muts <- c("ASXL1", "BCOR", "DNMT3A", "EZH2", "IDH2", "NRAS", "PHF6", "RUNX1", 
              "SRSF2", "U2AF1", "ZRSR2")

mut_hr_clust <- lapply(test_muts, function(gene){
  lapply(groups, function(cl){
    df <- subset(clinical_blasts, sub_group == cl)
    df_est <- getEstimates(gene, df) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_iws <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, clinical_mds) %>%
      mutate(sub_group = "IPSSM cohort")
  }) %>% Reduce(f = rbind)

mut_hr_comb <- rbind(mut_hr_clust, mut_hr_iws) %>%
  mutate(sub_group = factor(sub_group, levels = c(groups, "IPSSM cohort"))) %>%
  filter(N >= 10) %>%
  mutate(sub_group = droplevels(sub_group)) 

colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#0000FF", "#FFFFFF")

png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations.png", res = 300, heigh = 1300, width = 2200)
mut_hr_comb %>%
  ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-Group") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = colors) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())     
dev.off()

summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ delY*ASXL1 + AGE, clinical_mds))
summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ plus8*IDH2 + AGE, clinical_mds))


plotPair <- function(group, gene){
    df <- clinical_blasts %>%
        filter(sub_group == group) %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))
    
    p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ gene, df) %>%
        ggsurvplot(data = df, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 legend.labs  = c("mut", "wt")) +
        xlab("Time (Years)") 

    df2 <- clinical_mds %>%
        mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
        filter(!is.na(gene))
    p2 <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ gene, df2) %>%
        ggsurvplot(data = df2, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 legend.labs  = c("mut", "wt")) +
        xlab("Time (Years)")

    plot_grid(
        plot_grid(p$plot + 
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
            ylab("OS probability") +
            xlab("Time (years)") +
            ggtitle(paste(group, "vs", gene)) +
            theme(plot.title = element_text(hjust = 0.5)), 
            p$table, ncol = 1),
        plot_grid(p2$plot + 
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
            ylab("OS probability") +
            xlab("Time (years)") +
            ggtitle(paste("IPSSM cohort vs", gene)) +
            theme(plot.title = element_text(hjust = 0.5)), 
            p2$table, ncol = 1),
    ncol = 2, rel_widths = c(0.6, 1))
}

png("figures/GESMD_IWS_clustering/subgroup_interaction/TETm_BCOR.png", res = 300, heigh = 1500, width = 2800)
plotPair("TET2 monoallelic", "BCOR")
dev.off()


png("figures/GESMD_IWS_clustering/subgroup_interaction/TETm_EZH2.png", res = 300, heigh = 1500, width = 2800)
plotPair("TET2 monoallelic", "EZH2")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/TETbi_ASXL1.png", res = 300, heigh = 1500, width = 2800)
plotPair("TET2 bi-allelic", "ASXL1")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/TETbi_RUNX1.png", res = 300, heigh = 1500, width = 2800)
plotPair("TET2 bi-allelic", "RUNX1")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/plus8_RUNX1.png", res = 300, heigh = 1500, width = 2800)
plotPair("8+", "RUNX1")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/plus8_IDH2.png", res = 300, heigh = 1500, width = 2800)
plotPair("8+", "IDH2")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_RUNX1.png", res = 300, heigh = 1500, width = 2800)
plotPair("7-", "RUNX1")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/del20q_DNMT3A.png", res = 300, heigh = 1500, width = 2800)
plotPair("del20q", "DNMT3A")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_NRAS.png", res = 300, heigh = 1500, width = 2800)
plotPair("STAG2", "NRAS")
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_RUNX1.png", res = 300, heigh = 1500, width = 2800)
plotPair("STAG2", "RUNX1")
dev.off()