#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define prognosis for STAG2 mutated patients
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for subgroups exploration
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.7 R
#'
#' ---------------------------



# Load libraries and data
library(cowplot)
library(survminer)
library(survival)
library(tidyverse)
library(ipssm)
library(forestploter)

## Load data
load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/preprocess/clinical_preproc.Rdata")
load("data/GESMD/gesmd_data_1125.Rdata")
load("results/hershberger/hershberger_full.Rdata")
load("results/hershberger/hershberger_mds.Rdata")

classifySamples <- function(df){
    new_class <- case_when(
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2 bi-allelic",
    df$complex == 1     ~ "Complex",
    df$del7 == 1       ~ "7-",
    df$del5q == 1      ~ "del5q",
    df$SF3B1 == 1      ~ "SF3B1",
    df$STAG2 == 1       ~ "STAG2",
    df$BM_BLAST <= 5    ~ "Low blasts",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c("EZH2", "TET2 bi-allelic",  "7-", "STAG2", "del5q", "SF3B1", "Complex",
      "Low blasts", "MDS-IB1", "MDS-IB2"))
}

## Load data for full cohorts
IWS_full <- clinical %>%
    filter(!WHO_2016 %in% c("aCML", "CMML", "MDS/MPN-RS-T", "MDS/MPN-U", "other")) %>%
    mutate(PLT2 = pmin(PLT, 250))

gesmd <- gesmd_data_1125 %>%
    mutate(ID = as.character(register_number),
          WHO_2016 = who2017)

gesmd_full <- gesmd %>%
    filter(!WHO_2016 %in% c("WHO2017_LMA", "WHO2017_LMMC", "WHO2017_LMMC_0", 
        "WHO2017_LMMC_1", "WHO2017_LMMC_2", "WHO2017_LMMCX", "WHO2017_OTROS", 
        "WHO2017_SMD_SMP_INCLASIFICABLE", "WHO2017_SMD_SMP_SA_T", 
        "WHO2017_SMD_INCLASIFICABLE")) %>%
    filter(!is.na(TP53))

ipssm_process <- IPSSMprocess(gesmd_full)
ipssm_res <- IPSSMmain(ipssm_process)
ipssm_annot <- IPSSMannotate(ipssm_res)

gesmd_full <- left_join(gesmd_full, ipssm_annot %>% 
    mutate(IPSSM = IPSSMcat_mean,
    IPSSM = gsub(" ", "-", IPSSM , fixed = TRUE),
    IPSSM_SCORE = IPSSMscore_mean) %>%
    select(ID, IPSSM, IPSSM_SCORE), by = "ID")

joint_full <- bind_rows(
    IWS_full %>% mutate(dataset = "IWS") %>% mutate(complex = ifelse(complex == "complex", 1, 0)),
    gesmd_full %>% mutate(dataset = "GESMD") %>% filter(if_all(c(BM_BLAST, EZH2, STAG2, del7, TET2bi, del5q, SF3B1, complex), ~ !is.na(.)))) %>%
    mutate(sub_group = classifySamples(.),
    sub_group = ifelse(consensus %in% c("del5q", "mutated SF3B1", "Mutated TP53"), 
        consensus, as.character(sub_group)))
    
joint_mds <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS"),
    gesmd_dataset %>% mutate(dataset = "GESMD") 
)

joint_prognosis_plot <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS") ,
    gesmd_dataset %>% mutate(dataset = "GESMD")) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
     IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 


joint_prognosis <- joint_prognosis_plot %>% 
  mutate(sub_group = relevel(sub_group, ref = "Low blasts"), 
  mol_manual = ifelse(as.character(mol_manual) %in% c("Other", "CCUS-like"), "Low blasts", as.character(mol_manual)),
  mol_manual = relevel(factor(mol_manual), ref = "Low blasts"),
  IPSSM = factor(IPSSM, levels = c("Low", "Very-Low", "Moderate-Low", "Moderate-High", "High", "Very-High")))

hersh_prognosis <- hersh_mds %>%
  mutate(sub_group = relevel(sub_group, ref = "Low blasts"), 
  IPSSM = factor(IPSSM, levels = c("Low", "Very Low", "Moderate Low", "Moderate High", "High", "Very High")))


plotFitTable <- function(df, mod_var, title, palette, labels,ylab, xlim ){

    mod <- formula(paste("Surv(OS_YEARS, OS_STATUS) ~", mod_var))
    my_fit <- survfit(formula = mod, df) 
    my_fit$call$formula <- mod
    plot <- ggsurvplot(my_fit, data = df, surv.median.line = "hv", palette = palette,
        risk.table = TRUE, break.time.by = 2, conf.int = FALSE,
        xlim = xlim, ylim = c(0, 1),
        legend.labs  = labels)


    plot_grid(plot$plot + 
        ggtitle(title) +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab(ylab) +
         xlab("Time (years)"), 
        plot$table, ncol = 1)
}

plotPanelDatasets <- function(df, var_mod, title, palette, ylab = "OS prob.", xlim){

    labels <- levels(df[[var_mod]])
    IWS_df <- subset(df, dataset == "IWS")
    IWS_plot <- plotFitTable(IWS_df, var_mod, paste("IWS", title), palette, labels, ylab, xlim)

    GESMD_df <- subset(df, dataset == "GESMD")
    labels_gesmd <- labels[labels %in% unique(GESMD_df[[var_mod]])]
    palette_gesmd <- palette[labels %in% unique(GESMD_df[[var_mod]])]
    GESMD_plot <- plotFitTable(GESMD_df, var_mod, paste("GESMD", title), palette_gesmd, labels_gesmd, ylab, xlim)

    Hersh_df <- subset(df, dataset == "Hershberger")
    labels_hersh <- labels[labels %in% unique(Hersh_df[[var_mod]])] 
    palette_hersh <- palette[labels %in% unique(Hersh_df[[var_mod]])]
    Hersh_plot <- plotFitTable(Hersh_df, var_mod, paste("Hershberger", title), palette_hersh, labels_hersh, ylab, xlim)


    plot_grid(IWS_plot, GESMD_plot, Hersh_plot, ncol = 3)

}
createTitle <- function(title)  { 
    ggplot() + 
  labs(title = title) +
  theme_void() + 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
}
ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")
stag2_col <- c("#F5E4EE", "#EBCADB", "#D997BB", "#CC79A7", "#A15F83", "#663C53")

prognosis_stag2 <- bind_rows(
    joint_prognosis %>% filter(sub_group %in% c("STAG2", "Low blasts")),
    hersh_prognosis %>% filter(sub_group %in% c("STAG2", "Low blasts")) %>% 
    mutate(dataset = "Hershberger",
    IPSSM = fct_recode(IPSSM,
  "Very-Low" = "Very Low",
  "Moderate-Low" = "Moderate Low",
  "Moderate-High" = "Moderate High",
  "Very-High" = "Very High")
)) %>%
    mutate(sub_group = droplevels(sub_group),
        sub_group_short = fct_recode(sub_group,   "LB" = "Low blasts"),
        IPSSM_short = fct_recode(IPSSM, 
        "VL" = "Very-Low", "L" = "Low", "ML" = "Moderate-Low", 
        "MH" = "Moderate-High", "H" = "High", "VH" = "Very-High"),
    interaction = paste(sub_group_short, IPSSM_short, sep = " - "),
    interaction = factor(interaction, levels = c(paste(rep(c("LB", "STAG2"), each = 6), 
        c("VL", "L", "ML", "MH", "H", "VH"), sep = " - "))),
        dataset = factor(dataset, levels = c("IWS", "GESMD", "Hershberger"))) 


# stag2_ipssm <- plotPanelDatasets(prognosis_stag2, "interaction", "", c(ipssm_cols, stag2_col), xlim = c(0, 15))

stag2_ipssm <- plotPanelDatasets(filter(prognosis_stag2, sub_group =="STAG2"), "IPSSM_short", "", ipssm_cols, xlim = c(0, 15))

png("figures/GESMD_IWS_clustering/stag2/OS_interaction.png", width = 4500, height = 2000, res = 300)
plot_grid(createTitle("STAG2 - IPSSM interaction"), stag2_ipssm, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

prognosis_stag2b <- prognosis_stag2  %>%
    mutate(interaction = fct_recode(interaction,
    "STAG2 - VL-L-ML" = "STAG2 - VL", 
    "STAG2 - VL-L-ML" = "STAG2 - L", 
    "STAG2 - VL-L-ML" = "STAG2 - ML"))


stag2_ipssmb <- plotPanelDatasets(prognosis_stag2b, "interaction", "", 
    c(ipssm_cols, stag2_col[c(1, 3, 5, 6)]), xlim = c(0, 15))


png("figures/GESMD_IWS_clustering/stag2/OS_interaction_low.png", width = 4000, height = 2000, res = 300)
plot_grid(createTitle("STAG2 - IPSSM interaction"), stag2_ipssmb, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

prognosis_stag2_only <- prognosis_stag2 %>%
    filter(sub_group == "STAG2") %>%
    mutate(RUNX1 = ifelse(RUNX1 == 1, "Mutated RUNX1", "WT RUNX1"),
    RUNX1 = factor(RUNX1, levels = c("WT RUNX1", "Mutated RUNX1")),
     plus8 = ifelse(plus8 == 1, "8+", "WT"), 
     plus8 = factor(plus8, levels = c("WT","8+")),
     U2AF1 = ifelse(U2AF1 == 1, "Mutated U2AF1", "WT U2AF1"),
     U2AF1 = factor(U2AF1, levels = c("WT U2AF1", "Mutated U2AF1")),
     BM_BLAST_cat = case_when( BM_BLAST <= 5 ~ "<5%",
      BM_BLAST > 5 & BM_BLAST <= 10 ~ "5-10%", 
      BM_BLAST > 10 ~ "10+%", TRUE ~ NA_character_), 
      BM_BLAST_cat = factor(BM_BLAST_cat, levels = c("<5%", "5-10%", "10+%")))
      
stag2_runx1 <- plotPanelDatasets(prognosis_stag2_only, "RUNX1", "", c("#0072B2", "#D55E00"),  xlim = c(0, 15))
png("figures/GESMD_IWS_clustering/stag2/OS_RUNX1.png", width = 4000, height = 1500, res = 300)
plot_grid(createTitle("STAG2 - RUNX1"), stag2_runx1, ncol = 1, rel_heights = c(0.1, 1))
dev.off()



stag2_plus8 <- plotPanelDatasets(prognosis_stag2_only, "plus8", "", c("#0072B2", "#D55E00"), xlim = c(0, 15))
png("figures/GESMD_IWS_clustering/stag2/OS_plus8.png", width = 4000, height = 1500, res = 300)
plot_grid(createTitle("STAG2 - 8+"), stag2_plus8, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

stag2_u2af1 <- plotPanelDatasets(prognosis_stag2_only, "U2AF1", "", c("#0072B2", "#D55E00"), xlim = c(0, 15))
png("figures/GESMD_IWS_clustering/stag2/OS_u2af1.png", width = 4000, height = 1500, res = 300)
plot_grid(createTitle("STAG2 - U2AF1"), stag2_u2af1, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

stag2_bm <- plotPanelDatasets(prognosis_stag2_only, "BM_BLAST_cat", "", c("#D3D3D3", "#808080", "#1A1A1A"), xlim = c(0, 15))
png("figures/GESMD_IWS_clustering/stag2/OS_BM_BLASTS.png", width = 4000, height = 1500, res = 300)
plot_grid(createTitle("STAG2 - BM Blast %"), stag2_bm, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

### Interaction analysis
mod_joint_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ sub_group * IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
iws_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ sub_group * IPSSM_SCORE + AGE + SEX, data = filter(joint_prognosis, dataset == "IWS"))
gesmd_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ sub_group * IPSSM_SCORE + AGE + SEX, data = filter(joint_prognosis, dataset == "GESMD"))
hersh_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ sub_group * IPSSM_SCORE + AGE + SEX, data = hersh_prognosis)



stag2_mod_coefs <- lapply(list(iws_int, gesmd_int, mod_joint_int, hersh_int), function(mod) { 
    
    sum <- summary(mod)
    coefs <- sum$conf.int["sub_groupSTAG2:IPSSM_SCORE",]
}) %>%
Reduce(rbind,.)

n_samples_lb <- sapply(c("IWS", "GESMD"), function(x) 
    sum(joint_prognosis$dataset == x & joint_prognosis$sub_group == "Low blasts"))
n_samples_lb <- c(n_samples_lb, sum(n_samples_lb), sum(hersh_prognosis$sub_group == "Low blasts"))

n_samples_stag2 <- sapply(c("IWS", "GESMD"), function(x) 
    sum(joint_prognosis$dataset == x & joint_prognosis$sub_group == "STAG2"))
n_samples_stag2 <- c(n_samples_stag2, sum(n_samples_stag2), sum(hersh_prognosis$sub_group == "STAG2"))

df_for <- data.frame(Cohort = c("   IWS", "   GESMD", " Joint", "   Hershberger"), 
                        LB = n_samples_lb,
                        STAG2 = n_samples_stag2,
                        est = stag2_mod_coefs[, 1],
                        low = stag2_mod_coefs[, 3],
                        high = stag2_mod_coefs[, 4],
                        CI = paste(rep(" ", 40), collapse = "")
  )
df_for$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)",
                                         df_for$est, df_for$low, df_for$high)

df_for <- rbind(list("Discovery", "", "", NA, NA, NA, "", ""), df_for)
df_for <- rbind(list("Replication", "", "", NA, NA, NA, "", ""), df_for)

df_for <- df_for[c(1, 3:5, 2, 6), ]

ipsmm_forest <- forest(df_for[, c(1, 2, 3, 7, 8)],
              est = df_for$est,
              lower = df_for$low, 
              upper = df_for$high,
              sizes = 1,
              ci_column = 4,
              x_trans = "log" ,
              is_summary = c(rep(FALSE, 3), TRUE, FALSE, FALSE),
              title = "STAG2 - IPSSM interaction", 
              xlab = "Hazard Ratio", 
              theme = forest_theme(title_just = "center")) 
png("figures/GESMD_IWS_clustering/stag2/interaction_forest.png", width = 2000, height = 900, res = 300)
print(ipsmm_forest)
dev.off()


### IPSSM 2
n_samples_full <- sapply(list(IWS_full, gesmd_full, hersh_all_mds), nrow)
joint_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = prognosis_stag2_only)
iws_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "IWS"))
gesmd_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "GESMD"))
hersh_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "Hershberger"))

iwsfull_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_full)
gesmdfull_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = gesmd_full)
hershfull_ipssm <- coxph(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = hersh_all_mds)



ipssm_mod_coefs <- lapply(list(iws_ipssm, gesmd_ipssm, joint_ipssm, hersh_ipssm, iwsfull_ipssm, gesmdfull_ipssm, hershfull_ipssm), function(mod) { 
    
    sum <- summary(mod)
    coefs <- sum$conf.int["IPSSM_SCORE",]
}) %>%
Reduce(rbind,.)

n_samples_ipssm <- sapply(c("IWS", "GESMD", "Hershberger"), function(x) 
    sum(prognosis_stag2_only$dataset == x))
n_samples_ipssm <- c(n_samples_ipssm, n_samples_ipssm["IWS"] + n_samples_ipssm["GESMD"], n_samples_full)


df_ipssm <- data.frame(Cohort = c("   IWS", "   GESMD", " Joint", "   Hershberger", "   IWS", "   GESMD","   Hershberger"), 
                        N = n_samples_ipssm[c(1, 2, 4, 3, 5:7)],
                        est = ipssm_mod_coefs[, 1],
                        low = ipssm_mod_coefs[, 3],
                        high = ipssm_mod_coefs[, 4],
                        CI = paste(rep(" ", 40), collapse = "")
  )
df_ipssm$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)",
                                         df_ipssm$est, df_ipssm$low, df_ipssm$high)

df_ipssm <- rbind(list("Discovery", "",  NA, NA, NA, "", ""), df_ipssm)
df_ipssm <- rbind(list("Replication", "", NA, NA, NA, "", ""), df_ipssm)
df_ipssm <- rbind(list("Full cohort", "", NA, NA, NA, "", ""), df_ipssm)

df_ipssm <- df_ipssm[c(3:6, 2, 7, 1, 8:10), ]

ipsmm_forest2 <- forest(df_ipssm[, c(1, 2, 6, 7)],
              est = df_ipssm$est,
              lower = df_ipssm$low, 
              upper = df_ipssm$high,
              sizes = 1,
              ci_column = 3,
              x_trans = "log" ,
              is_summary = c(rep(FALSE, 3), TRUE, rep(FALSE, 6)),
              title = "STAG2 - IPSSM interaction", 
              xlab = "Hazard Ratio", 
              theme = forest_theme(title_just = "center")) 
png("figures/GESMD_IWS_clustering/stag2/interaction_forest2.png", width = 2000, height = 900, res = 300)
print(ipsmm_forest2)
dev.off()





## RUNX1
runx1_joint_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX + dataset, data = filter(prognosis_stag2_only, dataset != "Hershberger"))
iws_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "IWS"))
gesmd_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "GESMD"))
hersh_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "Hershberger"))

iwsfull_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = IWS_full)
gesmdfull_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = gesmd_full)
hershfull_runx1 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ RUNX1 + AGE + SEX, data = hersh_all_mds)

runx1_mod_coefs <- rbind(lapply(list(iws_runx1, gesmd_runx1, runx1_joint_int, hersh_runx1), function(mod) { 
    
    sum <- summary(mod)
    coefs <- sum$conf.int["RUNX1Mutated RUNX1",]
    }) %>%
    Reduce(rbind,.),
    lapply(list(iwsfull_runx1, gesmdfull_runx1, hershfull_runx1), function(mod) { 
    sum <- summary(mod)
    coefs <- sum$conf.int["RUNX1",]
    }) %>%
    Reduce(rbind,.)
)



n_samples_wtrunx1 <- sapply(c("IWS", "GESMD", "Hershberger"), function(x) 
    sum(prognosis_stag2_only$dataset == x & prognosis_stag2_only$RUNX1 == "WT RUNX1"))
n_samples_wtrunx1 <- c(n_samples_wtrunx1,n_samples_wtrunx1["IWS"] + n_samples_wtrunx1["GESMD"], 
    sum(IWS_full$RUNX1 == 0), sum(gesmd_full$RUNX1 == 0), sum(hersh_all_mds$RUNX1 == 0, na.rm = TRUE))

n_samples_mutrunx1 <- sapply(c("IWS", "GESMD", "Hershberger"), function(x) 
    sum(prognosis_stag2_only$dataset == x & prognosis_stag2_only$RUNX1 == "Mutated RUNX1"))
n_samples_mutrunx1 <- c(n_samples_mutrunx1,n_samples_mutrunx1["IWS"] + n_samples_mutrunx1["GESMD"],
    sum(IWS_full$RUNX1 == 1), sum(gesmd_full$RUNX1 == 1), sum(hersh_all_mds$RUNX1 == 1, na.rm = TRUE))




df_for_runx1 <- data.frame(Cohort = c("   IWS", "   GESMD", " Joint", "   Hershberger", "   IWS", "   GESMD","   Hershberger"), 
                        WT = n_samples_wtrunx1[c(1, 2, 4, 3, 5:7)],
                        RUNX1 = n_samples_mutrunx1[c(1, 2, 4, 3, 5:7)],
                        est = runx1_mod_coefs[, 1],
                        low = runx1_mod_coefs[, 3],
                        high = runx1_mod_coefs[, 4],
                        CI = paste(rep(" ", 40), collapse = "")
  )
df_for_runx1$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)",
                                         df_for_runx1$est, df_for_runx1$low, df_for_runx1$high)

df_for_runx1 <- rbind(list("Discovery", "", "", NA, NA, NA, "", ""), df_for_runx1)
df_for_runx1 <- rbind(list("Replication", "", "", NA, NA, NA, "", ""), df_for_runx1)
df_for_runx1 <- rbind(list("Full cohort", "", "", NA, NA, NA, "", ""), df_for_runx1)

df_for_runx1 <- df_for_runx1[c(3:6, 2, 7, 1, 8:10), ]

runx1_forest <- forest(df_for_runx1[, c(1, 2, 3, 7, 8)],
              est = df_for_runx1$est,
              lower = df_for_runx1$low, 
              upper = df_for_runx1$high,
              sizes = 1,
              ci_column = 4,
              x_trans = "log" ,
              is_summary = c(rep(FALSE, 3), TRUE, rep(FALSE, 6)),
              title = "STAG2 - RUNX1 effect", 
            ticks_at = c(0.5, 1, 2, 4),
              xlab = "Hazard Ratio", 
              theme = forest_theme(title_just = "center")) 
png("figures/GESMD_IWS_clustering/stag2/RUNX1_forest.png", width = 2000, height = 900, res = 300)
print(runx1_forest)
dev.off()

## 8+
plus8_joint_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX + dataset, data = filter(prognosis_stag2_only, dataset != "Hershberger"))
iws_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "IWS"))
gesmd_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "GESMD"))
hersh_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "Hershberger"))

iwsfull_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = IWS_full)
gesmdfull_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = gesmd_full)
hershfull_plus8 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ plus8 + AGE + SEX, data = hersh_all_mds)



plus8_mod_coefs <- rbind(lapply(list(iws_plus8, gesmd_plus8, plus8_joint_int, hersh_plus8), function(mod) { 
    
    sum <- summary(mod)
    coefs <- sum$conf.int["plus88+",]
    }) %>%
    Reduce(rbind,.),
    lapply(list(iwsfull_plus8, gesmdfull_plus8, hershfull_plus8), function(mod) { 
    sum <- summary(mod)
    coefs <- sum$conf.int["plus8",]
    }) %>%
    Reduce(rbind,.)
)


n_samples_wtplus8 <- sapply(c("IWS", "GESMD", "Hershberger"), function(x) 
    sum(prognosis_stag2_only$dataset == x & prognosis_stag2_only$plus8 == "WT"))
n_samples_wtplus8 <- c(n_samples_wtplus8,n_samples_wtplus8["IWS"] + n_samples_wtplus8["GESMD"], 
    sum(IWS_full$plus8 == 0), sum(gesmd_full$plus8 == 0, na.rm = TRUE), sum(hersh_all_mds$plus8 == 0, na.rm = TRUE))

n_samples_mutplus8 <- sapply(c("IWS", "GESMD", "Hershberger"), function(x) 
    sum(prognosis_stag2_only$dataset == x & prognosis_stag2_only$plus8 == "8+"))
n_samples_mutplus8 <- c(n_samples_mutplus8,n_samples_mutplus8["IWS"] + n_samples_mutplus8["GESMD"],
    sum(IWS_full$plus8 == 1), sum(gesmd_full$plus8 == 1, na.rm = TRUE), sum(hersh_all_mds$plus8 == 1, na.rm = TRUE))

df_for_plus8 <- data.frame(Cohort = c("   IWS", "   GESMD", " Joint", "   Hershberger", "   IWS", "   GESMD","   Hershberger"), 
                        WT = n_samples_wtplus8[c(1, 2, 4, 3, 5:7)],
                        plus8 = n_samples_mutplus8[c(1, 2, 4, 3, 5:7)],
                        est = plus8_mod_coefs[, 1],
                        low = plus8_mod_coefs[, 3],
                        high = plus8_mod_coefs[, 4],
                        CI = paste(rep(" ", 40), collapse = "")
  )
df_for_plus8$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)",
                                         df_for_plus8$est, df_for_plus8$low, df_for_plus8$high)

df_for_plus8 <- rbind(list("Discovery", "", "", NA, NA, NA, "", ""), df_for_plus8)
df_for_plus8 <- rbind(list("Replication", "", "", NA, NA, NA, "", ""), df_for_plus8)
df_for_plus8 <- rbind(list("Full cohort", "", "", NA, NA, NA, "", ""), df_for_plus8)

df_for_plus8 <- df_for_plus8[c(3:6, 2, 7, 1, 8:10), ]

plus8_forest <- forest(df_for_plus8[, c(1, 2, 3, 7, 8)],
              est = df_for_plus8$est,
              lower = df_for_plus8$low, 
              upper = df_for_plus8$high,
              sizes = 1,
              ci_column = 4,
              x_trans = "log" ,
              is_summary = c(rep(FALSE, 3), TRUE, rep(FALSE, 6)),
              title = "STAG2 - 8+ effect", 
            ticks_at = c(0.5, 1, 2, 4),
              xlab = "Hazard Ratio", 
              theme = forest_theme(title_just = "center")) 
png("figures/GESMD_IWS_clustering/stag2/plus8_forest.png", width = 2000, height = 900, res = 300)
print(plus8_forest)
dev.off()


## BM_BLASTS
bm_blasts_joint_int <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX + dataset, data = filter(prognosis_stag2_only, dataset != "Hershberger"))
iws_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "IWS"))
gesmd_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "GESMD"))
hersh_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = filter(prognosis_stag2_only, dataset == "Hershberger"))

iwsfull_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = IWS_full)
gesmdfull_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = gesmd_full)
hershfull_bm_blasts <- coxph(Surv(OS_YEARS, OS_STATUS) ~ BM_BLAST + AGE + SEX, data = hersh_all_mds)




bm_blasts_mod_coefs <- lapply(list(iws_bm_blasts, gesmd_bm_blasts, bm_blasts_joint_int, hersh_bm_blasts, iwsfull_bm_blasts, gesmdfull_bm_blasts, hershfull_bm_blasts), function(mod) { 
    
    sum <- summary(mod)
    coefs <- sum$conf.int["BM_BLAST",]
}) %>%
Reduce(rbind,.)

df_for_bm_blasts <- data.frame(Cohort = c("   IWS", "   GESMD", " Joint", "   Hershberger", "   IWS", "   GESMD","   Hershberger"), 
                        Samples = n_samples_ipssm[c(1, 2, 4, 3, 5:7)],
                        est = bm_blasts_mod_coefs[, 1],
                        low = bm_blasts_mod_coefs[, 3],
                        high = bm_blasts_mod_coefs[, 4],
                        CI = paste(rep(" ", 40), collapse = "")
  )
df_for_bm_blasts$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)",
                                         df_for_bm_blasts$est, df_for_bm_blasts$low, df_for_bm_blasts$high)

df_for_bm_blasts <- rbind(list("Discovery", "",  NA, NA, NA, "", ""), df_for_bm_blasts)
df_for_bm_blasts <- rbind(list("Replication", "", NA, NA, NA, "", ""), df_for_bm_blasts)
df_for_bm_blasts <- rbind(list("Full cohort", "", NA, NA, NA, "", ""), df_for_bm_blasts)

df_for_bm_blasts <- df_for_bm_blasts[c(3:6, 2, 7, 1, 8:10), ]

bm_blasts_forest <- forest(df_for_bm_blasts[, c(1, 2, 6, 7)],
              est = df_for_bm_blasts$est,
              lower = df_for_bm_blasts$low, 
              upper = df_for_bm_blasts$high,
              sizes = 1,
              ci_column = 3,
              x_trans = "log" ,
              is_summary = c(rep(FALSE, 3), TRUE, rep(FALSE, 6)),
              title = "STAG2 - BM BLAST (%) effect", 
            ticks_at = c(0.9, 1, 1.1, 1.2),
              xlab = "Hazard Ratio", 
              theme = forest_theme(title_just = "center")) 
png("figures/GESMD_IWS_clustering/stag2/bm_blasts_forest.png", width = 2000, height = 900, res = 300)
print(bm_blasts_forest)
dev.off()

png("figures/GESMD_IWS_clustering/stag2/stag2_panel.png", width = 5000, height = 4500, res = 300) 
plot_grid(stag2_ipssm, ipsmm_forest2,
    stag2_runx1, runx1_forest,
    stag2_plus8, plus8_forest,
    stag2_bm, bm_blasts_forest,
     ncol = 2, rel_widths = c(2, 1), rel_heights = c(1.2, 1, 1, 1),
     labels = "AUTO") 
dev.off()