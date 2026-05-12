#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compare prognosis between the subgroups of MDS
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for subgroups exploration
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
library(ipssm)
library(forestploter)


load("results/GESMD_IWS_clustering/gesmd_IWS_mds.Rdata")
load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/hershberger/hershberger_mds.Rdata")
load("results/hershberger/hershberger_full.Rdata")


## Remove 3 classical molecular groups
IWS_mds_f <- IWS_mds %>% 
  filter(!sub_group %in% c("Complex", "SF3B1-IB", "del5q-IB")) %>%
  mutate(sub_group = droplevels(sub_group))
gesmd_dataset_f <- gesmd_dataset %>% filter(!sub_group %in% c("Complex", "SF3B1-IB", "del5q-IB")) %>%
  mutate(sub_group = droplevels(sub_group))
hersh_mds_f <- hersh_mds %>% filter(!sub_group %in% c("Complex", "SF3B1-IB", "del5q-IB")) %>%
  mutate(sub_group = droplevels(sub_group))



joint_full <- bind_rows(
    IWS_full %>% mutate(dataset = "IWS") %>% mutate(complex = ifelse(complex == "complex", 1, 0)),
    gesmd_full %>% mutate(dataset = "GESMD") %>%
    mutate(CYTO_IPSSR = case_when(
        CYTO_IPSSR == "Intermediate" ~ "Int",
        CYTO_IPSSR == "Very Poor" ~ "Very-Poor",
        CYTO_IPSSR == "Very Good" ~ "Very-Good",
        TRUE ~ CYTO_IPSSR
    ),
    CYTO_IPSSR = factor(CYTO_IPSSR, levels = c( "Good", "Very-Good", "Int", "Poor", "Very-Poor"))
)
)
joint_full$CYTO_IPSSR = factor(joint_full$CYTO_IPSSR, levels = c( "Very-Good", "Good",  "Int", "Poor", "Very-Poor"))


joint_mds <- bind_rows(
    IWS_mds_f %>% mutate(dataset = "IWS"),
    gesmd_dataset_f %>% mutate(dataset = "GESMD") %>%
    mutate(CYTO_IPSSR = case_when(
        CYTO_IPSSR == "Intermediate" ~ "Int",
        CYTO_IPSSR == "Very Poor" ~ "Very-Poor",
        CYTO_IPSSR == "Very Good" ~ "Very-Good",
        TRUE ~ CYTO_IPSSR
    ),
    CYTO_IPSSR = factor(CYTO_IPSSR, levels = c( "Good", "Very-Good", "Int", "Poor", "Very-Poor"))
)
)
joint_mds$CYTO_IPSSR = factor(joint_mds$CYTO_IPSSR, levels = c( "Very-Good", "Good",  "Int", "Poor", "Very-Poor"))
## Overall survival
colors_all <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
    "#D55E00", "#999999", "grey40",  "black")
colors <-  c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
     "#999999", "grey40",  "black")

surv_IWS <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, IWS_mds_f) %>%
    ggsurvplot(data = IWS_mds_f, surv.median.line = "hv", palette = colors,
     risk.table = TRUE, break.time.by = 2, xlim = c(0, 10),
      legend.labs  = levels(IWS_mds_f$sub_group))

surv_gesmd <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ sub_group, gesmd_dataset_f) %>%
    ggsurvplot(data = gesmd_dataset_f, surv.median.line = "hv", 
               palette = colors, risk.table = TRUE, break.time.by = 2, 
      xlim = c(0, 10), legend.labs  = levels(gesmd_dataset_f$sub_group)) 

surv_hersh <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ sub_group, hersh_mds_f) %>%
    ggsurvplot(data = hersh_mds_f, surv.median.line = "hv", 
               palette = colors, risk.table = TRUE, break.time.by = 2, 
      xlim = c(0, 10), legend.labs  = levels(hersh_mds_f$sub_group))


os_comb <- plot_grid(
    plot_grid(surv_IWS$plot + 
        ggtitle("IWS") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("OS probability") +
         xlab("Time (years)"), 
        surv_IWS$table, ncol = 1),
    plot_grid(surv_gesmd$plot +
        ggtitle("GESMD") +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
        ylab("OS probability") +
         xlab("Time (years)"), 
        surv_gesmd$table, ncol = 1),
    plot_grid(surv_hersh$plot +
        ggtitle("MLL") +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
        ylab("OS probability") +
         xlab("Time (years)"), 
        surv_hersh$table, ncol = 1),
    ncol = 3)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_all_subgroups.png", width = 4500, height = 2000, res = 300)
os_comb
dev.off()

joint_prognosis_plot <- bind_rows(
    IWS_mds_f %>% mutate(dataset = "IWS") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), STAG2, IPSSM, IPSSM_SCORE, mol_manual, sub_group, AGE, SEX, dataset, BM_BLAST),
    gesmd_dataset_f %>% mutate(dataset = "GESMD") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), STAG2, IPSSM, IPSSM_SCORE,  mol_manual, sub_group, AGE, SEX, dataset, BM_BLAST)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
     IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 


joint_prognosis <- joint_prognosis_plot %>% 
  #mutate(sub_group = relevel(sub_group, ref = "MDS-LB"), 
  mutate(IPSSM = factor(IPSSM, levels = c("Low", "Very-Low", "Moderate-Low", "Moderate-High", "High", "Very-High")))


# raw <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group, data = joint_prognosis))
# mod1 <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_prognosis))
# main_os_sum <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) )
# int_os_sum <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) )
# summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + BM_BLAST + AGE + SEX + dataset, data = joint_prognosis) )

# int_os_sum2 <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, 
#   subset = !IPSSM %in% c("Very-Low", "Low")) )


# new_groups <- c("EZH2", "STAG2", "-7", "TET2-bi")
# int_models <- lapply(new_groups, function(group){
#   df <- joint_prognosis %>%
#     mutate(main_group = ifelse(sub_group == group, group, "Rest"),
#     main_group = factor(main_group, levels = c(group, "Rest")))
#     summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ main_group*IPSSM_SCORE + BM_BLAST + AGE + SEX + dataset, data = df) )

# })
# names(int_models) <- new_groups

# int_models2 <- lapply(new_groups, function(group){
#   df <- joint_prognosis %>%
#     mutate(main_group = ifelse(sub_group == group, group, "Rest"),
#     main_group = factor(main_group, levels = c(group, "Rest"))) %>%
#     filter(main_group == group | !sub_group %in% new_groups)
#     summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ main_group*IPSSM_SCORE + BM_BLAST + AGE + SEX + dataset, data = df) )
# })
# names(int_models2) <- new_groups


# makeModelTab <- function(cox_sum){

#   coefs <- cox_sum$coefficients
#   conf_int <- cox_sum$conf.int

#   tab <- tibble(
#     Variable = gsub("sub_group", "", rownames(coefs)),
#     HR = sprintf("%.1f (%.1f - %.1f)", 
#       round(coefs[, "exp(coef)"], 2), round(conf_int[, "lower .95"], 2), round(conf_int[, "upper .95"], 2)),
#     p_value = signif(coefs[, "Pr(>|z|)"], 2)
#   )
#   tab
# }

# mods <- list(raw = raw, mod1 = mod1, main_os_sum = main_os_sum, int_os_sum = int_os_sum)
# os_res <- lapply(names(mods), function(x) makeModelTab(mods[[x]]) %>% mutate(model = x)) %>% Reduce(f = rbind)


# os_res_filt <- os_res %>% 
#   filter(str_detect(Variable, paste(groups, collapse = "|")))
# write.table(os_res_filt, 
#             file = "results/GESMD_IWS_clustering/OS_subgroups_coxph.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             col.names = TRUE,
#             row.names = FALSE)
groups <- levels(joint_prognosis$sub_group)
names(groups) <- groups

main_model_tabs <- lapply(groups, function(group){

  joint_prognosis_sub <- joint_prognosis %>%
    filter(sub_group == group) %>%
    mutate(sub_group = droplevels(sub_group),
    IPSSM = relevel(IPSSM, ref = "Moderate-High"))
    mod <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX + dataset, data = joint_prognosis_sub))
  
  cat_sum <- table(joint_prognosis_sub$IPSSM)
  sel_cats <- cat_sum[cat_sum >= 10] %>% names()
  coefs <- mod$coefficients
  conf_int <- mod$conf.int
  tab <- tibble(
    Variable = gsub("IPSSM", "", rownames(coefs)),
    Group = group,
    HR = round(coefs[, "exp(coef)"], 2),
    HR_Low = round(conf_int[, "lower .95"], 2), 
    HR_High = round(conf_int[, "upper .95"], 2),
    p_value = signif(coefs[, "Pr(>|z|)"], 2)
  ) %>%
  filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
  filter(Variable %in% sel_cats)
  tab
}) %>% Reduce(f = rbind)


hr_os_group_plot <- bind_rows(main_model_tabs,
  tibble(Variable = "Moderate-High", Group = groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
ggplot(aes(x = Variable, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors) +
  geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "IPSSM",
    color = "Sub-group") +
  ggtitle("OS (Ref: IPSSM Moderate-High)")

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_os_group_plot
dev.off()


## Plot median OS
ipssm_raw_surv <- survfit( Surv(OS_YEARS,OS_STATUS) ~ sub_group+IPSSM, joint_prognosis)
ipssm_raw_surv_tib <- data.frame(summary(ipssm_raw_surv)$table) %>%
  rownames_to_column("Group") %>%
  separate(Group, into = c("sub_group", "IPSSM"), sep = "\\, ") %>%
  mutate(sub_group = gsub("sub_group=", "", sub_group),
  IPSSM = gsub("IPSSM=", "", IPSSM),
  IPSSM = gsub(" ", "", IPSSM)) %>%
  as_tibble() %>%
  mutate(sub_group = factor(sub_group, levels = groups),
  IPSSM = factor(IPSSM, c("Very-High", "High", "Moderate-High", "Moderate-Low", "Low", "Very-Low"))) %>%
  filter(records >= 10 & !is.na(median) ) %>%
  mutate(UCL = ifelse(is.na(X0.95UCL), median + 3 , X0.95UCL),)


median_OS_plot <- ggplot(ipssm_raw_surv_tib, aes(x = IPSSM, y = median, color = sub_group)) +
  geom_point(size = 2) +
  geom_segment(aes(y = `X0.95LCL`, yend = UCL),
  arrow = arrow(length = unit(ifelse(is.na(ipssm_raw_surv_tib$`X0.95UCL`), 0.3, 0), "cm"))) +
  coord_flip() +
  scale_color_manual(values = colors) +
  labs(title = "Median OS by sub-group and IPSSM",
       x = "", color = "",
       y = "Median OS (years)") +
  facet_grid(sub_group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_median.png", width = 1500, height = 2000, res = 300)
median_OS_plot
dev.off()


## Compute table for IPSSM effect for sub-group
getCoefs <- function(df, group, IPSSM_filter = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")){
  df_sub <- df %>%
    filter(sub_group == group & IPSSM %in% IPSSM_filter)
  mod <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = df_sub))
  mod$coefficients[1, c(1, 5)]
}
ipssm_effect_list <- lapply(groups, function(group) getCoefs(joint_prognosis, group))
ipssm_effect_list_lowrisk <- lapply(groups, function(group) 
  getCoefs(joint_prognosis, group, IPSSM_filter = c("Very-Low", "Low", "Moderate-Low")))
ipssm_effect_list_highrisk <- lapply(groups, function(group) 
  getCoefs(joint_prognosis, group, IPSSM_filter = c("Very-High", "High", "Moderate-High")))

ipssm_effect_list_joint <- c(ipssm_effect_list, ipssm_effect_list_lowrisk, ipssm_effect_list_highrisk)

ipssm_effect_tib <- tibble(group = names(ipssm_effect_list_joint), 
  IPSSM_effect = sapply(ipssm_effect_list_joint, function(x) x[1]),
  p_value = sapply(ipssm_effect_list_joint, function(x) x[2]),
  Prog = "OS", Category = rep(c("All", "Low-risk", "High-risk"), each = length(groups))) 

# joint_prognosis_os_plot <- joint_prognosis %>%
#     mutate(comb_group = paste(sub_group, IPSSM, sep = "_"),
#     comb_group = factor(comb_group, levels = unique(comb_group)),
#     comb_group = relevel(comb_group, ref = "Low blasts_Low"))

# cat_count <- table(joint_prognosis_os_plot$comb_group)
# sel_cats <- names(cat_count[cat_count >= 10])

# joint_prognosis_os_plot <- joint_prognosis_os_plot %>%
#   filter(comb_group %in% sel_cats)  %>%
#   mutate(comb_group = droplevels(comb_group))
# os_mod_subgroup_ipssm <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ comb_group + AGE + SEX + dataset, data = joint_prognosis_os_plot))

# os_mod_subgroup_ipssm_tab <- tibble(
#     Variable = gsub("comb_group", "", rownames(os_mod_subgroup_ipssm$coefficients)),
#     Group = sapply(strsplit(Variable, "_", 2), `[`, 1),
#     IPSSM = sapply(strsplit(Variable, "_", 2), `[`, 2),
#     HR = round(os_mod_subgroup_ipssm$coefficients[, "exp(coef)"], 2),
#     HR_Low = round(os_mod_subgroup_ipssm$conf.int[, "lower .95"], 2), 
#     HR_High = round(os_mod_subgroup_ipssm$conf.int[, "upper .95"], 2),
#     p_value = signif(os_mod_subgroup_ipssm$coefficients[, "Pr(>|z|)"], 2)
#   ) %>%
#   filter(IPSSM != "NA") %>%
#   filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
#   bind_rows(., tibble(
#     Variable = "Low blasts_Low",
#     Group = "Low blasts",
#     IPSSM = "Low",
#     HR = 1,
#     HR_Low = 1,
#     HR_High = 1,
#     p_value = NA
#   ))



# hr_os_group_plot2 <- os_mod_subgroup_ipssm_tab %>%
# mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
#       Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
# ggplot(aes(x = Group, y = HR, color = Group)) +
#   geom_point() +
#   scale_color_manual(values = colors6) +
#   geom_errorbar(aes(x = Group, ymin = HR_Low, ymax = HR_High)) +
#   theme_bw() +
#   facet_grid(. ~ IPSSM, scale = "free_x", space = "free_x") +
#   scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#   plot.title = element_text(hjust = 0.5)) +
#   labs(x = "Sub-group",
#   color = "Sub-group") +
#   ggtitle("OS (Ref: Low blasts/IPSSM Low)")

# png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
# hr_os_group_plot2
# dev.off()



## OS by IPSSM (IWS + GESMD)
IPSSM_groups <- levels(joint_prognosis_plot$IPSSM)
names(IPSSM_groups) <- IPSSM_groups

os_joint_IPSSM <- lapply(IPSSM_groups, function(cat){
  df <- filter(joint_prognosis_plot, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors[which(groups %in% sel_clusts)],
               legend.labs = groups[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2, xlim = c(0, 10)) +
    xlab("Time (Years)")
  p
})

os_ipssm_plots <- lapply(IPSSM_groups, function(ipssm){

    os_plot <-   plot_grid(
            plot_grid(os_joint_IPSSM[[ipssm]]$plot + 
                ggtitle(paste(ipssm, "IPSSM")) +
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                ylab("OS probability") +
                xlab("Time (years)"), 
                os_joint_IPSSM[[ipssm]]$table, ncol = 1),
            ncol = 1)

              
    ggsave(plot = os_plot, filename = paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_", ipssm, "_joint.png"),
        width = 1800, height = 2000, dpi = 300, units = "px")
    os_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_IPSSM_joint_panel.png", width = 4000, height = 6000, res = 300)
plot_grid(plotlist = os_ipssm_plots, ncol = 2, labels = "AUTO")
dev.off()



## Survival by subgroup 
ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")
survs_subgroups <- lapply(groups, function(group){
  df <- filter(joint_prognosis_plot, sub_group == group)
  n_ipssm <- df %>% 
    group_by(IPSSM) %>%
    summarize(n = n())
  sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM)
  df <- filter(df, IPSSM %in% sel_ipssm)
  p <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ IPSSM, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = ipssm_cols[which(IPSSM_groups %in% sel_ipssm)],
               legend.labs = IPSSM_groups[which(IPSSM_groups %in% sel_ipssm)],
               risk.table = TRUE, break.time.by = 2, xlim = c(0, 10)) +
    xlab("Time (Years)")
  p
})
names(survs_subgroups) <- groups

survs_subgroups_plots <- lapply(names(survs_subgroups), function(group){
    
   os_plot <- plot_grid(
        plot_grid(survs_subgroups[[group]]$plot + 
            ggtitle(group) +
             theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
             ylab("OS probability") +
             xlab("Time (years)"), 
            survs_subgroups[[group]]$table, ncol = 1),
        ncol = 1)
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_", group, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
    os_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_panel.png", width = 4000, height = 6000, res = 300)
plot_grid(plotlist = survs_subgroups_plots, ncol = 2, labels = "AUTO")
dev.off()

## AML transformation
IWS_full_amlt <- IWS_full %>%
    mutate(PROG_STATE = ifelse(AMLt_STATUS == 1, "AMLt", 
      ifelse(OS_STATUS == 1, "Death", "Censored")),
      PROG_STATE = factor(PROG_STATE, levels = c("Censored", "AMLt", "Death")))


amlt_IWS <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group, IWS_mds_f) %>%
    ggsurvplot(data = IWS_mds_f,  palette = colors,
     risk.table = TRUE, break.time.by = 2, fun = "event", xlim = c(0, 10),
      legend.labs  = levels(IWS_mds_f$sub_group))

aml_all <- plot_grid(amlt_IWS$plot + 
        ggtitle("IWS") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("AMLt probability") +
         xlab("Time (years)"), 
        amlt_IWS$table, ncol = 1)
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_all_subgroups.png", width = 2000, height = 2000, res = 300)
aml_all
dev.off()

# IWS_aml <- mutate(IWS_mds, sub_group = relevel(sub_group, ref = "Low blasts"), IPSSM = relevel(IPSSM, ref = "Low"))

# aml_raw <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group, data = IWS_aml))
# aml_adj <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group + AGE + SEX, data = IWS_aml))
# aml_main <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = IWS_aml) )
# aml_int <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX , data = IWS_aml) )


# aml_mods <- list(raw = aml_raw, adj = aml_adj, main_aml = aml_main, int_aml = aml_int)
# aml_res <- lapply(names(aml_mods), function(x) makeModelTab(aml_mods[[x]]) %>% mutate(model = x)) %>% Reduce(f = rbind)

# aml_res_filt <- aml_res %>% 
#   filter(str_detect(Variable, paste(groups, collapse = "|")))
# write.table(aml_res_filt, 
#             file = "results/GESMD_IWS_clustering/AMLt_subgroups_coxph.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             col.names = TRUE,
#             row.names = FALSE)

aml_model_tabs <- lapply(groups, function(group){

  IWS_sub <- IWS_mds_f %>%
    filter(sub_group == group) %>%
    mutate(sub_group = droplevels(sub_group),
    IPSSM = relevel(IPSSM, ref = "Moderate-High"))
    mod <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM + AGE + SEX, data = IWS_sub))
  
  cat_sum <- table(IWS_sub$IPSSM)
  sel_cats <- cat_sum[cat_sum >= 10] %>% names()
  coefs <- mod$coefficients
  conf_int <- mod$conf.int
  tab <- tibble(
    Variable = gsub("IPSSM", "", rownames(coefs)),
    Group = group,
    HR = round(coefs[, "exp(coef)"], 2),
    HR_Low = round(conf_int[, "lower .95"], 2), 
    HR_High = round(conf_int[, "upper .95"], 2),
    p_value = signif(coefs[, "Pr(>|z|)"], 2)
  ) %>%
  filter(!Variable %in% c("AGE", "SEXM")) %>%
  filter(Variable %in% sel_cats)
  tab
}) %>% Reduce(f = rbind)


hr_aml_group_plot <- bind_rows(aml_model_tabs,
  tibble(Variable = "Moderate-High", Group = groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
      filter(HR_Low != 0) %>%
      filter(Group != "EZH2") %>%
ggplot(aes(x = Variable, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors[-1]) +
  geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "IPSSM",
    color = "Sub-group") +
  ggtitle("AMLt (Reference: IPSSM Moderate-High)")

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_aml_group_plot
dev.off()

## Plot median AMLt
IWS_amlt <- IWS_mds_f %>%
    mutate(PROG_STATE = ifelse(AMLt_STATUS == 1, "AMLt", 
      ifelse(OS_STATUS == 1, "Death", "Censored")),
      PROG_STATE = factor(PROG_STATE, levels = c("Censored", "AMLt", "Death")))

ipssm_raw_amlt <- survfit( Surv(AMLt_YEARS,PROG_STATE) ~ sub_group + IPSSM, IWS_amlt, conf.type = "log-log")
ipssm_raw_amlt_sum <- summary(ipssm_raw_amlt, times = c(1, 2, 5), extend = TRUE)

lower <- ipssm_raw_amlt_sum$lower
colnames(lower) <- paste0("CILOW_", c("censored", "AMLt", "Death"))

upper <- ipssm_raw_amlt_sum$upper
colnames(upper) <- paste0("CIHIGH_", c("censored", "AMLt", "Death"))

pstate <- ipssm_raw_amlt_sum$pstate
colnames(pstate) <- paste0("P_", c("censored", "AMLt", "Death"))



ipssm_full_surv <- survfit( Surv(AMLt_YEARS,PROG_STATE) ~ IPSSM, IWS_amlt, conf.type = "log-log")
ipssm_full_amlt_sum <- summary(ipssm_full_surv, times = c(1, 2, 5), extend = TRUE)

lower_full <- ipssm_full_amlt_sum$lower
colnames(lower_full) <- paste0("CILOW_", c("censored", "AMLt", "Death"))

upper_full <- ipssm_full_amlt_sum$upper
colnames(upper_full) <- paste0("CIHIGH_", c("censored", "AMLt", "Death"))

pstate_full <- ipssm_full_amlt_sum$pstate
colnames(pstate_full) <- paste0("P_", c("censored", "AMLt", "Death"))



ipssm_raw_amlt_tib <- tibble(
  Group = ipssm_raw_amlt_sum$strata,
  Time = ipssm_raw_amlt_sum$time) %>%
  separate(Group, into = c("sub_group", "IPSSM"), sep = "\\, ") %>%
  mutate(sub_group = gsub("sub_group=", "", sub_group),
  IPSSM = gsub("IPSSM=", "", IPSSM),
  IPSSM = gsub(" ", "", IPSSM))  %>%
  bind_cols(., lower, upper, pstate) %>%
  pivot_longer(cols = starts_with(c("P_", "CI")), names_to = "var", values_to = "value") %>%
  separate(var, into = c("type", "Event"), sep = "_") %>%
  filter(Event == "AMLt") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  bind_rows(., 
  tibble( sub_group = "IWS full",
    IPSSM = ipssm_full_amlt_sum$strata,
  Time = ipssm_full_amlt_sum$time) %>%
  mutate(IPSSM = gsub("IPSSM=", "", IPSSM),
  IPSSM = gsub(" ", "", IPSSM))  %>%
  bind_cols(., lower_full, upper_full, pstate_full) %>%
  pivot_longer(cols = starts_with(c("P_", "CI")), names_to = "var", values_to = "value") %>%
  separate(var, into = c("type", "Event"), sep = "_") %>%
  filter(Event == "AMLt") %>%
  pivot_wider(names_from = type, values_from = value)) %>%
  mutate(time = factor(Time, levels = c("1", "2", "5")),
  IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
  sub_group = factor(sub_group, levels = c(groups, "IWS full")))

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_P_subgroups.png", width = 2500, height = 1200, res = 300)
# ipssm_raw_amlt_tib %>% 
#   filter(!is.na(P)) %>%
#   filter((sub_group == "EZH2" & IPSSM %in% c("Moderate-High", "High", "Very-High")) | 
#     sub_group == "TET2-bi" | sub_group == "MDS-LB" | (sub_group == "-7" & IPSSM %in% c("High", "Very-High")) | 
#     (sub_group == "STAG2" & IPSSM %in% c("Moderate-High", "High", "Very-High")) |
#     (sub_group == "MDS-IB1" & IPSSM != "Very-Low") |
#     (sub_group == "MDS-IB2" & IPSSM %in% c("Moderate-Low", "Moderate-High", "High", "Very-High"))) %>%
#   ggplot(aes(x = sub_group, y = P, fill = sub_group, color = sub_group)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymin = CILOW, ymax = CIHIGH)) +
#   facet_grid(time ~ IPSSM, scales = "free_x", space = "free_x") +
#   scale_fill_manual(values = c(colors, "black")) +
#   scale_color_manual(values =  c(colors, "black")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#   plot.title = element_text(hjust = 0.5)) +
#   labs(x = "Sub-group",
#     fill = "Sub-group", color = "Sub-group") +
#   ggtitle("AMLt probability at 1, 2, and 5 years")
# dev.off()


AMLt_P_plot <- ipssm_raw_amlt_tib %>% 
mutate(IPSSM = factor(IPSSM, levels = rev(levels(IPSSM)))) %>%
  filter(!is.na(P) & Time != 5) %>%
  filter((sub_group == "EZH2" & IPSSM %in% c("Moderate-High", "High", "Very-High")) | 
    sub_group == "TET2-bi" | sub_group == "MDS-LB" | (sub_group == "-7" & IPSSM %in% c("High", "Very-High")) | 
    (sub_group == "STAG2" & IPSSM %in% c("Moderate-High", "High", "Very-High")) |
    (sub_group == "MDS-IB1" & IPSSM != "Very-Low") |
    (sub_group == "MDS-IB2" & IPSSM %in% c("Moderate-Low", "Moderate-High", "High", "Very-High"))) %>%
    mutate(Time = factor(Time, levels = c("1", "2"), labels = c("1 year", "2 years"))) %>%
  ggplot(aes(x = IPSSM, y = P, fill = sub_group, color = sub_group)) +
  geom_point(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = CILOW, ymax = CIHIGH)) +
  facet_grid(sub_group ~ Time, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(colors, "black")) +
  scale_color_manual(values =  c(colors, "black")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sub-group", y = "Probability of AMLt",
    fill = "", color = "") +
  ggtitle("AMLt probability")

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_P_subgroups2.png", width = 2000, height = 2500, res = 300)
AMLt_P_plot
dev.off()


## Compute effect of AML
coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX,  data = IWS_full_amlt)

getCoefsAML <- function(df, group, IPSSM_cats = levels(df$IPSSM)){
  df_sub <- df %>%
    filter(sub_group == group & IPSSM %in% IPSSM_cats) 
  mod <- summary(coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX, data = df_sub))
  mod$coefficients[1, c(1, 5)]
}
ipssm_amlt_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt)
ipssm_amlt_lowRisk_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt, 
  IPSSM = c("Very-Low", "Low", "Moderate-Low"))
ipssm_amlt_HighRisk_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt, 
  IPSSM = c("Moderate-High", "High", "Very-High"))

ipssm_amlt_effect_list_joint <- c(ipssm_amlt_effect_list, ipssm_amlt_lowRisk_effect_list, ipssm_amlt_HighRisk_effect_list)

ipssm_amlt_effect_tib <- tibble(group = names(ipssm_amlt_effect_list_joint), 
  IPSSM_effect = sapply(ipssm_amlt_effect_list_joint, function(x) x[1]),
  p_value = sapply(ipssm_amlt_effect_list_joint, function(x) x[2]),
  Prog = "AMLt", Category = rep(c("All", "Low-risk", "High-risk"), each = length(groups)))

ipssm_prog_joint <- bind_rows(ipssm_effect_tib, ipssm_amlt_effect_tib)
write.table(ipssm_prog_joint, file = "results/GESMD_IWS_clustering/IPSSM_effects_OS_AMLt.txt", 
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX,  
  data = IWS_full_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))
coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX,  
  data = IWS_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))

IWS_full_amlt$TET2bi_group <- ifelse(IWS_full_amlt$ID %in% filter(IWS_amlt, sub_group == "TET2-bi")$ID, "TET2bi", "Other")
IWS_amlt$TET2bi_group <- ifelse(IWS_amlt$TET2bi == 1, "TET2bi", "Other")

IWS_amlt <- IWS_amlt %>%
  mutate(AMLt_STATUS2 = ifelse(AMLt_YEARS > 5, 0, AMLt_STATUS),
        AMLt_YEARS2 = ifelse(AMLt_YEARS > 5, 5, AMLt_YEARS),
      PROG_STATE2 = ifelse(AMLt_STATUS2 == 1, "AMLt", 
      ifelse(OS_STATUS == 1 & OS_YEARS < 5, "Death", "Censored")),
      PROG_STATE2 = factor(PROG_STATE2, levels = c("Censored", "AMLt", "Death")))


IWS_full_amlt <- IWS_full_amlt %>%
  mutate(AMLt_STATUS2 = ifelse(AMLt_YEARS > 5, 0, AMLt_STATUS),
        AMLt_YEARS2 = ifelse(AMLt_YEARS > 5, 5, AMLt_YEARS),
      PROG_STATE2 = ifelse(AMLt_STATUS2 == 1, "AMLt", 
      ifelse(OS_STATUS == 1 & OS_YEARS < 5, "Death", "Censored")),
      PROG_STATE2 = factor(PROG_STATE2, levels = c("Censored", "AMLt", "Death")))


coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))

coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_amlt)

coxph(Surv(AMLt_YEARS2,PROG_STATE2 == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))

coxph(Surv(AMLt_YEARS2,PROG_STATE2 == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_amlt)


coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_full_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))

coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_full_amlt)

coxph(Surv(AMLt_YEARS2,PROG_STATE2 == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_full_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))

coxph(Surv(AMLt_YEARS2,PROG_STATE2 == "AMLt") ~ TET2bi_group + IPSSM_SCORE + AGE + SEX,  
  data = IWS_full_amlt)


## TET2-bi effect in low-risk IPSSM
amlt_TET2bi_morph <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ TET2bi_group, 
  IWS_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low")) %>%
    ggsurvplot(data = IWS_amlt,  palette = c("black", "#56B4E9"),
     risk.table = TRUE, break.time.by = 2, fun = "event", xlim = c(0, 6),
     ylim = c(0, 0.2),
      legend.labs  = levels(IWS_amlt$TET2bi_group))

aml_tet2bi_morph <- plot_grid(amlt_TET2bi_morph$plot + 
        ggtitle("IWS Morphologic Low Risk") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("AMLt probability") +
         xlab("Time (years)"), 
        amlt_TET2bi_morph$table, ncol = 1)
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_tet2bi_morph_lowrisk.png", width = 2000, height = 2000, res = 300)
aml_tet2bi_morph
dev.off()

amlt_TET2bi_full <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ TET2bi_group, 
  IWS_full_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low")) %>%
    ggsurvplot(data = IWS_full_amlt,  palette = c("black", "#56B4E9"),
     risk.table = TRUE, break.time.by = 2, fun = "event", xlim = c(0, 6),
     ylim = c(0, 0.2), legend.labs = c("Rest", "TET2-bi"),
     pval = sprintf("HR = %.2f\nP = %.3f", 2.530043, 0.00457),
    pval.coord = c(0.1, 0.15))

aml_tet2bi_full <- plot_grid(amlt_TET2bi_full$plot + 
        ggtitle("IWS Full Low Risk") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("AMLt probability") +
         xlab("Time (years)"), 
        amlt_TET2bi_full$table + xlab("Time (years)"), ncol = 1, rel_heights = c(1.6, 1))
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_tet2bi_full_lowrisk.png", width = 2000, height = 2000, res = 300)
aml_tet2bi_full
dev.off()


IWS_full_amlt2 <- mutate(IWS_full_amlt, Risk = ifelse(IPSSM %in% c("Very-Low", "Low", "Moderate-Low"), "Low-risk", 
  ifelse(IPSSM %in% c("Moderate-High", "High", "Very-High"), "High-risk", NA)))

amlt_TET2bi_full2 <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ TET2bi_group + Risk, 
  IWS_full_amlt2) %>%
    ggsurvplot(data = IWS_full_amlt2,  
     risk.table = TRUE, break.time.by = 2, fun = "event", xlim = c(0, 6),
     ylim = c(0, 0.2))


# IWS_aml_plot <- IWS_mds %>%
#     mutate(comb_group = paste(sub_group, IPSSM, sep = "_"),
#     comb_group = factor(comb_group, levels = unique(comb_group)),
#     comb_group = relevel(comb_group, ref = "Low blasts_Low"))

# cat_count_aml <- table(IWS_aml_plot$comb_group)
# sel_cats_aml <- names(cat_count_aml[cat_count_aml >= 10])
# sel_cats_aml <- sel_cats_aml[!sel_cats_aml %in% c("EZH2_Moderate-High", "Low blasts_NA", "TET2 bi-allelic_Very-Low")]
# IWS_aml_plot <- IWS_aml_plot %>%
#   filter(comb_group %in% sel_cats_aml)  %>%
#   mutate(comb_group = droplevels(comb_group))
# aml_mod_subgroup_ipssm <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ comb_group + AGE + SEX, data = IWS_aml_plot))

# aml_mod_subgroup_ipssm_tab <- tibble(
#     Variable = gsub("comb_group", "", rownames(aml_mod_subgroup_ipssm$coefficients)),
#     Group = sapply(strsplit(Variable, "_", 2), `[`, 1),
#     IPSSM = sapply(strsplit(Variable, "_", 2), `[`, 2),
#     HR = round(aml_mod_subgroup_ipssm$coefficients[, "exp(coef)"], 2),
#     HR_Low = round(aml_mod_subgroup_ipssm$conf.int[, "lower .95"], 2), 
#     HR_High = round(aml_mod_subgroup_ipssm$conf.int[, "upper .95"], 2),
#     p_value = signif(aml_mod_subgroup_ipssm$coefficients[, "Pr(>|z|)"], 2)
#   ) %>%
#   filter(IPSSM != "NA") %>%
#   filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
#   bind_rows(., tibble(
#     Variable = "Low blasts_Low",
#     Group = "Low blasts",
#     IPSSM = "Low",
#     HR = 1,
#     HR_Low = 1,
#     HR_High = 1,
#     p_value = NA
# #   ))


# hr_aml_group_plot2 <- aml_mod_subgroup_ipssm_tab %>%
# mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
#       Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
#       filter(!IPSSM == "Very-Low") %>%
# ggplot(aes(x = Group, y = HR, color = Group)) +
#   geom_point() +
#   scale_color_manual(values = colors6) +
#   geom_errorbar(aes(x = Group, ymin = HR_Low, ymax = HR_High)) +
#   theme_bw() +
#   facet_grid(. ~ IPSSM, scale = "free_x", space = "free_x") +
#   scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), 2)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#   plot.title = element_text(hjust = 0.5)) +
#   labs(x = "Sub-group",
#   color = "Sub-group") +
#   ggtitle("AMLt (Ref: Low blasts/IPSSM Low)")

# png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
# hr_aml_group_plot2
# dev.off()

# panel_prognosis <- 
# plot_grid(
#   plot_grid(
#     os_comb,
#     plot_grid(hr_os_group_plot2, hr_os_group_plot, ncol = 1, labels = c("B", "C")),
#     ncol = 2, labels = c("A", "")),
#   plot_grid(
#     aml_all,
#     plot_grid(hr_aml_group_plot2, hr_aml_group_plot, ncol = 1, labels = c("E", "F")),
#     ncol = 2, labels = c("D", "")
#   ),
#   ncol = 1, labels = "")

# png("figures/GESMD_IWS_clustering/subgroup_prognosis/prognosis_panel.png", width = 4500, height = 4000, res = 300)
# panel_prognosis
# dev.off()


## AMLt by IPSSM (only IWS)
amlt_ipssm <- lapply(IPSSM_groups, function(cat){
  df <- filter(IWS_mds_f, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
               palette = colors[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2, xlim = c(0, 10)) +
    xlab("Time (Years)")
  p
})
amlt_ipssm_plots <- lapply(IPSSM_groups, function(ipssm){
   aml_plot <-  plot_grid(
        plot_grid(amlt_ipssm[[ipssm]]$plot + 
            ggtitle(ipssm) +
             theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
             ylab("AMLt probability") +
             xlab("Time (years)") +
             ylim(0, 1), 
            amlt_ipssm[[ipssm]]$table, ncol = 1),
        ncol = 1)
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_", ipssm, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
    aml_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_IPSSM_panel.png", width = 4000, height = 6000, res = 300)
plot_grid(plotlist = amlt_ipssm_plots, ncol = 2, labels = "AUTO")
dev.off()
## AMlt by subgroup (only IWS)
amlt_subgroups <- lapply(groups, function(group){
  df <- filter(IWS_mds, sub_group == group)
  n_ipssm <- df %>% 
    group_by(IPSSM) %>%
    summarize(n = n())
  sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM)
  df <- filter(df, IPSSM %in% sel_ipssm)
  p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ IPSSM, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
               palette = ipssm_cols[which(IPSSM_groups %in% sel_ipssm)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})
names(amlt_subgroups) <- groups
amlt_subgroups_plots <- lapply(names(amlt_subgroups), function(group){
    aml_plot <- plot_grid(
        plot_grid(amlt_subgroups[[group]]$plot + 
            ggtitle(group) +
             theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
             ylab("AMLt probability") +
             xlab("Time (years)") +
             ylim(0, 1), 
            amlt_subgroups[[group]]$table, ncol = 1),
        ncol = 1)
    png_filename <- paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_", group, ".png")
    png_filename <- gsub("[^A-Za-z0-9/_-].", "_", png_filename, fixed = TRUE)
    ggsave(png_filename, width = 2000, height = 2000, dpi = 300, units = "px")
    aml_plot
})
png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_panel.png", width = 4000, height = 6000, res = 300)
plot_grid(plotlist = amlt_subgroups_plots, ncol = 2, labels = "AUTO")
dev.off()


## Interactions
#####################################################################

## Compute HR for mutations and sub-groups
getEstimates <- function(gene, df, outcome = "OS", joint = FALSE){
  if (outcome == "OS"){
    out = "OS_YEARS,OS_STATUS"
  } 
  if (outcome == "AML"){
    out = "AMLt_YEARS,PROG_STATE == \"AMLt\""
  }

  formula_base <- paste("Surv(", out, ") ~", gene, "+ AGE + SEX")
  if (joint){
    formula_base <- paste(formula_base, "+ dataset")
  }

  mod <- summary(coxph(formula(formula_base), df))
  coef <- mod$coefficients
  conf <- mod$conf.int
  hr <- conf[1, 1]
  hr_l <- conf[1, 3]
  hr_h <- conf[1, 4]
  pval <- coef[1, 5]
  freq <- mean(df[[gene]], na.rm = TRUE)
  N <- sum(df[[gene]], na.rm = TRUE)
  data.frame(HR = hr, HR_L = hr_l, HR_H = hr_h, Pvalue = pval, 
             Gene = gene, Freq = freq, N = N )
}

## Select genes frequent in new sub-groups
test_muts <- c("ASXL1", "SRSF2", "DNMT3A", "TP53mono", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "CBL", "IDH1", "PHF6", "plus8", "delY", "del20q",
    "BM_BLAST", "HB", "PLT")

# mut_hr_clust <- lapply(test_muts, function(gene){
#   lapply(groups, function(cl){
#     df <- subset(IWS_mds, sub_group == cl)
#     df_est <- getEstimates(gene, df) %>%
#       mutate(sub_group = cl)
#   }) %>% Reduce(f = rbind)
# }) %>% Reduce(f = rbind)
# mut_hr_full <- lapply(test_muts, function(gene){
#   df_est <- getEstimates(gene, IWS_full) %>%
#       mutate(sub_group = "IPSSM cohort")
#   }) %>% Reduce(f = rbind)

# mut_hr_comb <- rbind(mut_hr_clust, mut_hr_full) %>%
#   mutate(sub_group = factor(sub_group, levels = c(groups, "IPSSM cohort"))) %>%
#   filter(N >= 10) %>%
#   mutate(sub_group = droplevels(sub_group)) %>%
#   as_tibble() %>%
#   filter(!Gene %in% c("DDX41", "IDH1")) %>%
#   mutate(Gene = factor(case_when(
#     Gene == "plus8" ~ "8+",
#     Gene == "TET2other" ~ "TET2 mono-allelic",
#     TRUE ~ Gene
#   )))


# os_int_plot <- mut_hr_comb %>%
#   ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-group") +
#   ggtitle("Overall Survival in IWS") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   scale_color_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5))    
# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations.png", res = 300, heigh = 1300, width = 2200)
#  os_int_plot
# dev.off()


mut_hr_clust_joint <- lapply(test_muts, function(gene){
  lapply(groups, function(cl){
    df <- subset(joint_mds, sub_group == cl)
    df_est <- getEstimates(gene, df, joint = TRUE) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_full_joint <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, joint_full, joint = TRUE) %>%
      mutate(sub_group = "Both cohorts")
  }) %>% Reduce(f = rbind)


cyto_hr <- lapply(groups, function(cl){
    df <- subset(joint_mds, sub_group == cl)
    df_est <- getEstimates("as.numeric(CYTO_IPSSR)", df, joint = TRUE) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind) %>%
  bind_rows(., getEstimates("as.numeric(CYTO_IPSSR)", joint_full, joint = TRUE) %>%
      mutate(sub_group = "Both cohorts")
  ) %>%
  mutate(Gene = "CYTO_IPSSR", N = 100)




mut_hr_comb_joint <- rbind(mut_hr_clust_joint, mut_hr_full_joint, cyto_hr) %>%
  mutate(sub_group = factor(sub_group, levels = c(groups, "Both cohorts"))) %>%
  filter(N >= 10) %>%
  mutate(sub_group = droplevels(sub_group)) %>%
  as_tibble() %>%
  mutate(Gene = factor(case_when(
    Gene == "plus8" ~ "+8",
    Gene == "delY" ~ "-Y",
    Gene == "TP53mono" ~ "TP53 mono-allelic",
    Gene == "TET2other" ~ "TET2 mono-allelic",
    Gene == "CYTO_IPSSR" ~ "Cytogenetic Risk",
    TRUE ~ Gene
  )),
  Gene = factor(Gene, levels = c("RUNX1", "IDH2", "CBL", "U2AF1", "SRSF2",
    "DNMT3A", "ASXL1", "BCOR", "IDH1", "PHF6", "SETBP1",
    "TP53 mono-allelic", "TET2 mono-allelic", "ZRSR2", "Cytogenetic Risk", "+8", "-Y", 
    "del20q",  "BM_BLAST", "HB", "PLT")))


os_int_sig <- data.frame(
  Gene = c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk", "Cytogenetic Risk"), 
  x = c("EZH2", "STAG2", "EZH2", "STAG2", "EZH2", "STAG2"),           
  y = c(4),                      
  label = c("*", "*", "*", "*", "*", "*")  
)



os_int_plot_main <- mut_hr_comb_joint %>%
  filter(!is.na(Gene) & Gene %in% c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk")) %>%
  ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-group") +
  ggtitle("Overall Survival") +
  facet_wrap(~ Gene, scales = "free_x") +
  geom_text(data = os_int_sig, aes(x = x, y = y, label = label), 
    color = "black", inherit.aes = FALSE) +
  scale_fill_manual(name = "", values = c(colors, "#59758a")) +
  scale_color_manual(name = "", values = c(colors, "#59758a")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  plot.title = element_text(hjust = 0.5), 
  legend.position = "bottom")    
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/HR_mutations_joint.png", res = 300, heigh = 1300, width = 2200)
 os_int_plot_main
dev.off()

os_int_plot_sub <- mut_hr_comb_joint %>%
  filter(!is.na(Gene) & !Gene %in% c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk", "BM_BLAST", "HB", "PLT")) %>%
  ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-group") +
  ggtitle("Overall Survival") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = c(colors, "#59758a")) +
  scale_color_manual(name = "", values = c(colors, "#59758a")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  plot.title = element_text(hjust = 0.5))    
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/HR_mutations_joint_sup.png", res = 300, heigh = 1300, width = 2200)
 os_int_plot_sub
dev.off()

mut_hr_comb2 <- mut_hr_comb_joint %>%
 group_by(Gene) %>% 
 mutate(HR_Ref = HR[sub_group == "Both cohorts"],
  HR_L_Ref = HR_L[sub_group == "Both cohorts"],
        HR_H_Ref = HR_H[sub_group == "Both cohorts"]) %>%
 ungroup()
top_low <- filter(mut_hr_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
    arrange(HR_H - HR_L_Ref)
top_high <- filter(mut_hr_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
    arrange(HR_H_Ref - HR_L)
filter(mut_hr_comb2, HR_H < HR_L_Ref | HR_L > HR_H_Ref)

filter(mut_hr_comb2, HR_H < HR_Ref | HR_L > HR_Ref)


joint_full_subs <- mutate(joint_full, 
  ezh2_group = ifelse(ID %in% filter(joint_mds, sub_group == "EZH2")$ID, "EZH2", "Other"),
  ezh2_group = factor(ezh2_group, levels = c("EZH2", "Other")),
  stag2_group = ifelse(ID %in% filter(joint_mds, sub_group == "STAG2")$ID, "STAG2", "Other"),
  stag2_group = factor(stag2_group, levels = c("STAG2", "Other")),
  mono7 = ifelse(ID %in% filter(joint_mds, sub_group == "-7")$ID, "-7", "Other"),
  mono7 = factor(mono7, levels = c("-7", "Other"))

)
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ ezh2_group*TET2other + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ stag2_group*U2AF1 + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ ezh2_group*SETBP1 + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ stag2_group*plus8 + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ mono7*BM_BLAST + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ stag2_group*BM_BLAST + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ ezh2_group*PLT + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ ezh2_group*as.numeric(CYTO_IPSSR) + AGE + SEX + dataset, joint_full_subs))
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ stag2_group*as.numeric(CYTO_IPSSR) + AGE + SEX + dataset, joint_full_subs))

summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ as.numeric(CYTO_IPSSR) + AGE + SEX + dataset, subset = sub_group == "STAG2", data = joint_mds))

stag2_cyto <- 
    survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ CYTO_IPSSR, joint_mds, 
    subset = sub_group == "STAG2" & CYTO_IPSSR %in% c("Good", "Int")) %>%
        ggsurvplot(data = joint_mds, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 legend = "none", title = "Cytogenetic Risk in STAG2 subgroup",
                  pval = sprintf("HR = %.2f\nP = %.3f", 0.994970, 0.98201),
                 pval.coord = c(7, 0.9), legend.labs = c("Good", "Int"),
                 xlim = c(0, 10), palette = c("#66bd63", "#fee08b" )
                 ) +
        xlab("Time (Years)") 

stag2_cyto$plot <- stag2_cyto$plot +  theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_cyto.png", width = 2000, height = 2000, res = 300)
stag2_cyto
dev.off()

full_cyto <- 
    survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ CYTO_IPSSR, joint_full) %>%
        ggsurvplot(data = joint_full, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, legend = "none",
                 title = "Cytogenetic Risk in full cohorts",
                 legend.labs = c("Very Good", "Good", "Int", "Poor", "Very Poor"),
                 xlim = c(0, 10), palette = c("#2ca25f", "#66bd63", "#fee08b", "#f46d43", "#d73027")) +
        xlab("Time (Years)") 
full_cyto$plot <- full_cyto$plot +  theme(plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/full_cyto.png", width = 2000, height = 2000, res = 300)
full_cyto
dev.off()



### Interaction analysis
fit_os_models_by_cohort <- function(var_name, data_list, outcome = "OS") {

  if (outcome == "OS") {
    time_event <- "Surv(OS_YEARS, OS_STATUS)"
  } else if (outcome == "AMLt") {
    time_event <- "Surv(AMLt_YEARS, PROG_STATE == 'AMLt')"
  } else {
    stop("Invalid outcome specified. Use 'OS' or 'AMLt'.")
  }
  joint_df <- data_list$Joint
  hersh_df <- data_list$Hershberger
  formula_joint <- formula(paste(time_event, "~", var_name, "+ AGE + SEX + dataset"))
  formula_single <- formula(paste(time_event, "~", var_name, "+ AGE + SEX"))
  formula_nosex <- formula(paste(time_event, "~", var_name, "+ AGE"))

  if (sum(hersh_df$SEX == "M") > 0 & sum(hersh_df$SEX == "F") > 0) {
    hersh_formula <- formula_single
  } else {
    hersh_formula <- formula_nosex
  }

  list(
    IWS = coxph(formula_single, data = filter(joint_df, dataset == "IWS")),
    GESMD = coxph(formula_single, data = filter(joint_df, dataset == "GESMD")),
    Joint = coxph(formula_joint, data = joint_df),
    MLL = coxph(hersh_formula, data = hersh_df)
  )
}

extract_term <- function(mod, term = 1) {
  mod_sum <- summary(mod)
  tibble(
    term = rownames(mod_sum$conf.int)[term],
    est = mod_sum$conf.int[term, "exp(coef)"],
    low = mod_sum$conf.int[term, "lower .95"],
    high = mod_sum$conf.int[term, "upper .95"]
  )
}

build_forest_table <- function(models, data_list, col1_fun, col2_fun, headers, term = 1) {
  model_stats <- lapply(models, extract_term, term = term) %>% bind_rows()

  sample_sizes <- lapply(data_list, function(df) {
    tibble(
      col1 = col1_fun(df),
      col2 = col2_fun(df)
    )
  }) %>% bind_rows()
  colnames(sample_sizes) <- headers

  df_for <- tibble(
    Cohort = names(models),
    est = model_stats$est,
    low = model_stats$low,
    high = model_stats$high,
    CI = strrep(" ", 40)
  ) %>%
  bind_cols(sample_sizes)
  df_for <- df_for[, c("Cohort", headers,  "CI", "est", "low", "high")]

  df_for$`HR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", df_for$est, df_for$low, df_for$high)

    df_for <- bind_rows(
      tibble(Cohort = "Discovery"),
      df_for[1:3, ],
      tibble(Cohort = "Replication"),
      df_for[4:nrow(df_for), ]
    )

  h1 <- df_for[, headers[1], drop = TRUE]
  h1 <- as.character(h1)
  h1[is.na(h1)] <- ""
  df_for[, 2] <- h1

  h2 <- df_for[, headers[2], drop = TRUE]
  h2 <- as.character(h2)
  h2[is.na(h2)] <- ""
  df_for[, 3] <- h2

  df_for$`HR (95% CI)`[is.na(df_for$`HR (95% CI)`)] <- ""
  df_for$CI[is.na(df_for$CI)] <- ""

  df_for
}

make_forest_table <- function(var_name, data_list, outcome = "OS", col1_fun, col2_fun, headers) {
  models <- fit_os_models_by_cohort(var_name = var_name, data_list = data_list, outcome = outcome)
  build_forest_table(models = models, data_list = data_list, col1_fun = col1_fun, col2_fun = col2_fun, headers = headers)
}

plot_forest <- function(df_for, plot_title, summary_rows, ci_column = 4) {


  forest_obj <- forest(
    df_for[, !colnames(df_for) %in% c("est", "low", "high")],
    est = df_for$est,
    lower = df_for$low,
    upper = df_for$high,
    sizes = 1,
    ci_column = ci_column,
    x_trans = "log",
    is_summary = summary_rows,
    title = plot_title,
    xlab = "Hazard Ratio",
    theme = forest_theme(title_just = "center"), 
  )
}

hersh_mds_cyto <- mutate(hersh_mds, 
CYTO_IPSSR = case_when(
    CYTO_IPSSR == "Intermediate" ~ "Int",
    TRUE ~ gsub(" ", "-", CYTO_IPSSR)
),
CYTO_IPSSR = factor(CYTO_IPSSR, levels = c("Very-Good", "Good", "Int", "Poor", "Very-Poor")))


hersh_all_mds_cyto <- mutate(hersh_all_mds, 
CYTO_IPSSR = case_when(
    CYTO_IPSSR == "Intermediate" ~ "Int",
    TRUE ~ gsub(" ", "-", CYTO_IPSSR)
),
CYTO_IPSSR = factor(CYTO_IPSSR, levels = c("Very-Good", "Good", "Int", "Poor", "Very-Poor")))



data_list_stag2 <- list(
    IWS = filter(joint_mds, dataset == "IWS" & sub_group == "STAG2"),
    GESMD = filter(joint_mds, dataset == "GESMD" & sub_group == "STAG2"),
    Joint = filter(joint_mds, sub_group == "STAG2"),
    Hershberger = filter(hersh_mds_cyto, sub_group == "STAG2")
)
data_list_full <- list(
    IWS = filter(joint_full, dataset == "IWS"),
    GESMD = filter(joint_full, dataset == "GESMD"),
    Joint = joint_full,
    Hershberger = hersh_all_mds_cyto %>% mutate(TET2other = TET2mono)
)
stag2_plus8 <- make_forest_table(var_name = "plus8", data_list = data_list_stag2, 
  col1_fun = function(df) sum(df$plus8 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$plus8 == 1, na.rm = TRUE),
  headers = c("WT", "+8"))

summary_rows <- rep(FALSE, 6)
summary_rows[4] <- TRUE


stag2_cyto_tab <- make_forest_table(var_name = "as.numeric(CYTO_IPSSR)", data_list = data_list_stag2, 
  col1_fun = function(df) sum(!is.na(df$CYTO_IPSSR)), 
  col2_fun = function(df) return(NA),
  headers = c("N CYTO", "A"))

stag2_cyto_forest <- plot_forest(stag2_cyto_tab[, -3], plot_title = "Cytogenetic Risk in STAG2 sub-group",
   summary_rows = summary_rows, ci_column = 3)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/forest_cyto_stag2.png", width = 2000, height = 2000, res = 300)
print(stag2_cyto_forest)
dev.off()

full_cyto_tab <- make_forest_table(var_name = "as.numeric(CYTO_IPSSR)", data_list = data_list_full, 
  col1_fun = function(df) sum(!is.na(df$CYTO_IPSSR)), 
  col2_fun = function(df) return(NA),
  headers = c("N CYTO", "A"))

full_cyto_forest <- plot_forest(full_cyto_tab[, -3], plot_title = "Cytogenetic Risk in full dataset",
   summary_rows = summary_rows, ci_column = 3)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/forest_cyto_full.png", width = 2000, height = 2000, res = 300)
print(full_cyto_forest)
dev.off()



createSurvFit <- function(df, variable, legend.labs, palette){

    model <- as.formula(paste("Surv(OS_YEARS,OS_STATUS) ~", variable))
    surv <- surv_fit(formula = model, df)
    ggsurvplot(fit = surv, data = df, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 palette = palette,
                 legend.labs  = legend.labs, xlim = c(0, 10))
}
makeSurvPlot <- function(df, variable, legend.labs = c("WT", "MUT"), palette = c("black", "red"), title){
    ggsurv <- createSurvFit(df, variable, legend.labs = legend.labs, palette = palette)
    plot_grid(ggsurv$plot + 
                theme(legend.position = "none", 
                    plot.title = element_text(hjust = 0.5)) +
                ylab("OS probability") +
                xlab("Time (years)") +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5)), 
                ggsurv$table + xlab("Time (years)"), ncol = 1, rel_heights = c(1, 0.4))

}

data_list_ezh2 <- list(
    IWS = filter(joint_mds, dataset == "IWS" & sub_group == "EZH2"),
    GESMD = filter(joint_mds, dataset == "GESMD" & sub_group == "EZH2"),
    Joint = filter(joint_mds, sub_group == "EZH2"),
    Hershberger = filter(hersh_mds_cyto, sub_group == "EZH2") %>% mutate(TET2other = TET2mono)
)

## EZH2 vs TET2 mono-allelic
ezh2_tet2m_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("EZH2")), variable = "TET2other", title = "TET2 mono-allelic in EZH2 subgroup")
ezh2_tet2m_tab <- make_forest_table(var_name = "TET2other", data_list = data_list_ezh2, 
  col1_fun = function(df) sum(df$TET2other == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$TET2other == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

ezh2_tet2m_forest <- plot_forest(ezh2_tet2m_tab, plot_title = "TET2 mono-allelic in EZH2 sub-group",
   summary_rows = summary_rows)

full_tet2m_surv <- makeSurvPlot(joint_full, variable = "TET2other", title = "TET2 mono-allelic in full cohort")
full_tet2m_tab <- make_forest_table(var_name = "TET2other", data_list = data_list_full, 
  col1_fun = function(df) sum(df$TET2other == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$TET2other == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

full_tet2m_forest <- plot_forest(full_tet2m_tab, plot_title = "TET2 mono-allelic in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/EZH2_TET2mono_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(ezh2_tet2m_surv, full_tet2m_surv, ezh2_tet2m_forest, full_tet2m_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()

## STAG2 vs U2AF1
stag2_u2af1_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("STAG2")), variable = "U2AF1", title = "U2AF1 in STAG2 subgroup")
stag2_u2af1_tab <- make_forest_table(var_name = "U2AF1", data_list = data_list_stag2, 
  col1_fun = function(df) sum(df$U2AF1 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$U2AF1 == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

stag2_u2af1_tab[3, c("est", "low", "high")] <- NA

stag2_u2af1_forest <- plot_forest(stag2_u2af1_tab, plot_title = "U2AF1 in STAG2 sub-group",
   summary_rows = summary_rows)

full_u2af1_surv <- makeSurvPlot(joint_full, variable = "U2AF1", title = "U2AF1 in full cohort")
full_u2af1_tab <- make_forest_table(var_name = "U2AF1", data_list = data_list_full, 
  col1_fun = function(df) sum(df$U2AF1 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$U2AF1 == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

full_u2af1_forest <- plot_forest(full_u2af1_tab, plot_title = "U2AF1 in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_U2AF1_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(stag2_u2af1_surv, full_u2af1_surv, stag2_u2af1_forest, full_u2af1_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()

## EZH2 vs SETBP1
ezh2_setbp1_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("EZH2")), variable = "SETBP1", title = "SETBP1 in EZH2 subgroup")
ezh2_setbp1_tab <- make_forest_table(var_name = "SETBP1", data_list = data_list_ezh2, 
  col1_fun = function(df) sum(df$SETBP1 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$SETBP1 == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

ezh2_setbp1_forest <- plot_forest(ezh2_setbp1_tab[1:4,], plot_title = "SETBP1 in EZH2 sub-group",
   summary_rows = summary_rows[1:4])

full_setbp1_surv <- makeSurvPlot(joint_full, variable = "SETBP1", title = "SETBP1 in full cohort")
full_setbp1_tab <- make_forest_table(var_name = "SETBP1", data_list = data_list_full, 
  col1_fun = function(df) sum(df$SETBP1 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$SETBP1 == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

full_setbp1_forest <- plot_forest(full_setbp1_tab, plot_title = "SETBP1 in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/EZH2_SETBP1_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(ezh2_setbp1_surv, full_setbp1_surv, ezh2_setbp1_forest, full_setbp1_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()

## STAG2 vs +8
stag2_plus8_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("STAG2")), variable = "plus8", 
  title = "+8 in STAG2 subgroup", legend.labs = c("WT", "+8"))
stag2_plus8_tab <- make_forest_table(var_name = "plus8", data_list = data_list_stag2, 
  col1_fun = function(df) sum(df$plus8 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$plus8 == 1, na.rm = TRUE), 
  headers = c("WT", "+8"))


stag2_plus8_forest <- plot_forest(stag2_plus8_tab, plot_title = "+8 in STAG2 sub-group",
   summary_rows = summary_rows)

full_plus8_surv <- makeSurvPlot(joint_full, variable = "plus8", 
  title = "+8 in full cohort", legend.labs = c("WT", "+8"))
full_plus8_tab <- make_forest_table(var_name = "plus8", data_list = data_list_full, 
  col1_fun = function(df) sum(df$plus8 == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$plus8 == 1, na.rm = TRUE), 
  headers = c("WT", "+8"))

full_plus8_forest <- plot_forest(full_plus8_tab, plot_title = "+8 in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_plus8_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(stag2_plus8_surv, full_plus8_surv, stag2_plus8_forest, full_plus8_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()

## -7 & STAG2 vs BM_BLAST
colores_blasts <- c("#c9bb92", "#FFB300", "#D32F2F")
data_list_mono7 <- list(
    IWS = filter(joint_mds, dataset == "IWS" & sub_group == "-7"),
    GESMD = filter(joint_mds, dataset == "GESMD" & sub_group == "-7"),
    Joint = filter(joint_mds, sub_group == "-7"),
    Hershberger = filter(hersh_mds, sub_group == "-7") 
)

mono7_blast_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("-7")) %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), 
  variable = "BLAST", title = "BM BLAST in -7 subgroup", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts)
mono7_blast_tab <- make_forest_table(var_name = "BM_BLAST", data_list = data_list_mono7,
  col1_fun = function(df) sum(!is.na(df$BM_BLAST)), 
  col2_fun = function(df) NA, 
  headers = c("N BM BLAST", "A"))
mono7_blast_forest <- plot_forest(mono7_blast_tab[, -3], plot_title = "BM BLAST in -7 subgroup",
   summary_rows = summary_rows, ci_column = 3)


stag2_blast_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("STAG2")) %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), 
  variable = "BLAST", title = "BM BLAST in STAG2 subgroup", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts)
stag2_blast_tab <- make_forest_table(var_name = "BM_BLAST", data_list = data_list_stag2,
  col1_fun = function(df) sum(!is.na(df$BM_BLAST)), 
  col2_fun = function(df) NA, 
  headers = c("N BM BLAST", "A"))
stag2_blast_forest <- plot_forest(stag2_blast_tab[, -3], plot_title = "BM BLAST in STAG2 subgroup",
   summary_rows = summary_rows, ci_column = 3)


full_blast_surv <- makeSurvPlot(joint_full %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), variable = "BLAST", 
  title = "BM BLAST in full cohort", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts)
full_blast_tab <- make_forest_table(var_name = "BM_BLAST", data_list = data_list_full, 
  col1_fun = function(df) sum(!is.na(df$BM_BLAST)), 
  col2_fun = function(df) NA, 
  headers = c("N BM BLAST", "A"))

full_blast_forest <- plot_forest(full_blast_tab, plot_title = "BM BLAST in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_mono7_blast_panel.png", width = 3500, height = 5000, res = 300)
plot_grid(
  plot_grid(stag2_blast_surv, mono7_blast_surv, stag2_blast_forest, mono7_blast_forest, 
    ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D")),
  plot_grid(full_blast_surv, full_blast_forest, ncol = 1, rel_heights = c(2, 1), labels = c("E", "F")),
  ncol = 1, rel_heights = c(1, 1)
)
dev.off()


## EZH2 vs PLT
colores_plt <- c("#C8E6C9", "#81C784", "#43A047", "#1B5E20")
ezh2_plt_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("EZH2")) %>% mutate(plt = cut(PLT, breaks = c(0, 50, 100, 150, Inf))), 
  variable = "plt", title = "PLT in EZH2 subgroup", legend.labs = c("<50", "50-100", "100-150", "\\>150"), 
  palette = colores_plt)
ezh2_plt_tab <- make_forest_table(var_name = "pmin(PLT, 250)", data_list = data_list_ezh2, 
  col1_fun = function(df) sum(!is.na(df$PLT)), 
  col2_fun = function(df) NA, 
  headers = c("N PLT", "A"))

ezh2_plt_tab <- mutate(ezh2_plt_tab, 
  `HR (95% CI)` = ifelse(is.na(est), "", sprintf("%.3f (%.3f to %.3f)", est, low, high))
)

ezh2_plt_forest <- plot_forest(ezh2_plt_tab[, -3], plot_title = "PLT in EZH2 sub-group",
   summary_rows = summary_rows, ci_column = 3)

full_plt_surv <- makeSurvPlot(joint_full %>% mutate(plt = cut(PLT, breaks = c(0, 50, 100, 150, Inf))), 
  variable = "plt", title = "PLT in full cohort", legend.labs = c("<50", "50-100", "100-150", "\\>150"), 
  palette = colores_plt)


full_plt_tab <- make_forest_table(var_name = "pmin(PLT, 250)", data_list = data_list_full, 
  col1_fun = function(df) sum(!is.na(df$PLT)), 
  col2_fun = function(df) NA, 
  headers = c("N PLT", "A"))
full_plt_tab <- mutate(full_plt_tab, 
  `HR (95% CI)` = ifelse(is.na(est), "", sprintf("%.3f (%.3f to %.3f)", est, low, high))
)

full_plt_forest <- plot_forest(full_plt_tab[, -3], plot_title = "PLT in full cohort",
   summary_rows = summary_rows, ci_column = 3)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/EZH2_PLT_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(ezh2_plt_surv, full_plt_surv, ezh2_plt_forest, full_plt_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()


## EZH2 vs CYTO
ezh2_cyto_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("EZH2") & CYTO_IPSSR != "Very-Good"), 
  variable = "CYTO_IPSSR", title = "Cytogenetic Risk in EZH2 subgroup", legend.labs = c("Good", "Int", "Poor"), 
  palette = c("#66bd63", "#fee08b", "#f46d43"))
ezh2_cyto_tab <- make_forest_table(var_name = "as.numeric(CYTO_IPSSR)", data_list = data_list_ezh2, 
  col1_fun = function(df) sum(!is.na(df$CYTO_IPSSR)), 
  col2_fun = function(df) NA, 
  headers = c("N Cytogenetic", "A"))

ezh2_cyto_forest <- plot_forest(ezh2_cyto_tab[, -3], plot_title = "Cytogenetic Risk in EZH2 sub-group",
   summary_rows = summary_rows, ci_column = 3)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/EZH2_CYTO_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(ezh2_cyto_surv, plot_grid(full_cyto$plot, full_cyto$table, ncol = 1, rel_heights = c(1.6, 1)),
  ezh2_cyto_forest, full_cyto_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
dev.off()


## Interactions in AMLt
mut_hr_amlt_clust_joint <- lapply(test_muts, function(gene){
  lapply(groups, function(cl){
    df <- subset(IWS_amlt, sub_group == cl)
    df_est <- getEstimates(gene, df, joint = FALSE, outcome = "AML") %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_amlt_full_joint <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, IWS_amlt, joint = FALSE, outcome = "AML") %>%
      mutate(sub_group = "IWS full")
  }) %>% Reduce(f = rbind)


cyto_hr_aml <- lapply(groups, function(cl){
    df <- subset(IWS_amlt, sub_group == cl)
    df_est <- getEstimates("as.numeric(CYTO_IPSSR)", df, joint = FALSE, outcome = "AML") %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind) %>%
  bind_rows(., getEstimates("as.numeric(CYTO_IPSSR)", IWS_amlt, joint = FALSE, outcome = "AML") %>%

      mutate(sub_group = "IWS full")
  ) %>%
  mutate(Gene = "CYTO_IPSSR", N = 100)


mut_hr_aml_comb_joint <- rbind(mut_hr_amlt_clust_joint, mut_hr_amlt_full_joint, cyto_hr_aml) %>%
  mutate(sub_group = factor(sub_group, levels = c(groups, "IWS full"))) %>%
  filter(N >= 10) %>%
  mutate(sub_group = droplevels(sub_group)) %>%
  as_tibble() %>%
  mutate(Gene = factor(case_when(
    Gene == "plus8" ~ "+8",
    Gene == "delY" ~ "-Y",
    Gene == "TP53mono" ~ "TP53 mono-allelic",
    Gene == "TET2other" ~ "TET2 mono-allelic",
    Gene == "CYTO_IPSSR" ~ "Cytogenetic Risk",
    TRUE ~ Gene
  )),
  Gene = factor(Gene, levels = c("RUNX1", "IDH2", "CBL", "U2AF1", "SRSF2",
    "DNMT3A", "ASXL1", "BCOR", "IDH1", "PHF6", "SETBP1",
    "TP53 mono-allelic", "TET2 mono-allelic", "ZRSR2", "Cytogenetic Risk", "+8", "-Y", 
    "del20q",  "BM_BLAST", "HB", "PLT")))




mut_hr_aml_comb2 <- mut_hr_aml_comb_joint %>%
 group_by(Gene) %>% 
 mutate(HR_Ref = HR[sub_group == "IWS full"],
  HR_L_Ref = HR_L[sub_group == "IWS full"],
        HR_H_Ref = HR_H[sub_group == "IWS full"]) %>%
 ungroup()
filter(mut_hr_aml_comb2, HR_H < HR_Ref | HR_L > HR_Ref) %>%
  arrange(sub_group)

IWS_aml_groups <- mutate(IWS_amlt, ezh2_group = ifelse(sub_group == "EZH2", 1, 0),
  stag2_group = ifelse(sub_group == "STAG2", 1, 0),
  mono7_group = ifelse(sub_group == "-7", 1, 0))

summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ ezh2_group*ASXL1 + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ ezh2_group*TET2other + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ ezh2_group*RUNX1 + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ mono7_group*BM_BLAST + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*ASXL1 + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*SRSF2 + AGE + SEX, IWS_aml_groups)) ## No signif
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*RUNX1 + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*BCOR + AGE + SEX, IWS_aml_groups))
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*IDH2 + AGE + SEX, IWS_aml_groups)) # No signif


amlt_int_sig <- data.frame(
  Gene = c("ASXL1", "TET2 mono-allelic", "RUNX1", "ASXL1", "RUNX1", "BCOR"), 
  x = c("EZH2", "EZH2", "EZH2", "STAG2", "STAG2", "STAG2"),           
  y = c(0.15),                      
  label = c("*", "**", "**", "*", "***", "*")  
)

aml_int_plot_main <- mut_hr_aml_comb_joint %>%
  filter(!is.na(Gene) & Gene %in% c("TET2 mono-allelic", "ASXL1", "RUNX1", "BCOR")) %>%
  ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-group") +
  ggtitle("AMLt") +
  facet_wrap(~ Gene, scales = "free_x") +
  geom_text(data = amlt_int_sig, aes(x = x, y = y, label = label), 
    color = "black", inherit.aes = FALSE) +
  scale_fill_manual(name = "", values = c(colors, "#59758a")) +
  scale_color_manual(name = "", values = c(colors, "#59758a")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  plot.title = element_text(hjust = 0.5), 
  legend.position = "bottom")    
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/HR_mutations_aml_joint.png", res = 300, heigh = 1300, width = 2200)
 aml_int_plot_main
dev.off()

aml_int_plot_sup <- mut_hr_aml_comb_joint %>%
  filter(!is.na(Gene) & !Gene %in% c("TET2 mono-allelic", "ASXL1", "RUNX1", "BCOR", "-Y", "BM_BLAST", "HB", "PLT", "TP53 mono-allelic", "del20q")) %>%
  filter(!(Gene == "DNMT3A" & sub_group == "-7")) %>%
  ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-group") +
  ggtitle("AMLt") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = c(colors, "#59758a")) +
  scale_color_manual(name = "", values = c(colors, "#59758a")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  plot.title = element_text(hjust = 0.5), 
  legend.position = "bottom")    
png("figures/GESMD_IWS_clustering/subgroup_prog_inter/HR_mutations_aml_sup.png", res = 300, heigh = 1300, width = 2200)
 aml_int_plot_sup
dev.off()

createSurvFitAML <- function(df, variable, legend.labs, palette, pval_height){

    base_model <- paste("Surv(AMLt_YEARS, PROG_STATE == 'AMLt') ~", variable)
    cox_model <- paste(base_model, "+ AGE + SEX")

    coxph_res <- summary(coxph(as.formula(cox_model), df))
    hr <- coxph_res$coefficients[1, 2]
    pval <- coxph_res$coefficients[1, 5]

    surv <- surv_fit(formula = as.formula(base_model), df)
    ggsurvplot(fit = surv, data = df, 
                 risk.table = TRUE, break.time.by = 2, 
                 palette = palette, fun = "event",
                 pval = sprintf("HR = %.2f\nP = %.3f", hr, pval),
                 pval.coord = c(3.5, pval_height),
                 legend.labs  = legend.labs, xlim = c(0, 5), ylim = c(0, 0.25))
}


makeSurvPlotAML <- function(df, variable, legend.labs = c("WT", "MUT"), palette = c("black", "red"), title, pval_height = 0.2){
    ggsurv <- createSurvFitAML(df, variable, legend.labs = legend.labs, palette = palette, pval_height = pval_height)
    plot_grid(ggsurv$plot + 
                theme(legend.position = "none", 
                    plot.title = element_text(hjust = 0.5)) +
                ylab("AMLt probability") +
                xlab("Time (years)") +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5)), 
                ggsurv$table + xlab("Time (years)"), ncol = 1, rel_heights = c(1.6, 1))
}



png("figures/GESMD_IWS_clustering/subgroup_prog_inter/prognosis_panel.png", width = 3200, height = 5500, res = 300)
plot_grid(
  plot_grid(median_OS_plot, AMLt_P_plot, 
  ncol = 1, rel_heights = c(1, 1.2), labels = c("A", "E")),
  plot_grid(os_int_plot_main,
      plot_grid(stag2_cyto$plot, stag2_cyto$table, ncol = 1, rel_heights = c(1.6, 1)),
  stag2_cyto_forest, aml_tet2bi_full, aml_int_plot_main, ncol = 1,
  rel_heights = c(1.3, 1.5, 1, 1.5, 1.3), labels = c("B", "C","D", "F", "G")),
  ncol = 2, rel_widths = c(1, 1.2)
)
dev.off()

## EZH2 vs TET2 mono
ezh2_tet2mono_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "EZH2"), 
  variable = "TET2other", title = "TET2 mono-allelic in EZH2 sub-group", pval_height = 0.10)

full_tet2mono_aml <- makeSurvPlotAML(IWS_amlt, 
  variable = "TET2other", title = "TET2 mono-allelic in IWS cohort", pval_height = 0.10)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/EZH2_TET2mono_aml_panel.png", res = 300, heigh = 1300, width = 2600)
plot_grid(ezh2_tet2mono_aml, full_tet2mono_aml, ncol = 2, labels = "AUTO")
dev.off()

## EZH2 and STAG2 vs ASXL1
ezh2_asxl1_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "EZH2"), 
  variable = "ASXL1", title = "ASXL1 in EZH2 sub-group", pval_height = 0.10)

stag2_asxl1_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "STAG2"), 
  variable = "ASXL1", title = "ASXL1 in STAG2 sub-group", pval_height = 0.10)

full_asxl1_aml <- makeSurvPlotAML(IWS_amlt, 
  variable = "ASXL1", title = "ASXL1 in IWS cohort", pval_height = 0.10)


png("figures/GESMD_IWS_clustering/subgroup_prog_inter/ASXL1_aml_panel.png", res = 300, height = 2000, width = 2600)
plot_grid(
    plot_grid(ezh2_asxl1_aml, stag2_asxl1_aml, nrow = 1, labels = c("A", "B")),
     full_asxl1_aml, ncol = 1, labels = c("", "C"))
dev.off()

## STAG2 and BCOR
stag2_bcor_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "STAG2"), 
  variable = "BCOR", title = "BCOR in STAG2 sub-group", pval_height = 0.10)

full_bcor_aml <- makeSurvPlotAML(IWS_amlt, 
  variable = "BCOR", title = "BCOR in IWS cohort", pval_height = 0.10)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_BCOR_aml_panel.png", res = 300, heigh = 1300, width = 2600)
plot_grid(stag2_bcor_aml, full_bcor_aml, ncol = 2, labels = "AUTO")
dev.off()

## EZH2 and STAG2 in RUNX1
ezh2_runx1_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "EZH2"), 
  variable = "RUNX1", title = "RUNX1 in EZH2 sub-group", pval_height = 0.10)

stag2_runx1_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "STAG2"), 
  variable = "RUNX1", title = "RUNX1 in STAG2 sub-group", pval_height = 0.10)

full_runx1_aml <- makeSurvPlotAML(IWS_amlt, 
  variable = "RUNX1", title = "RUNX1 in IWS cohort", pval_height = 0.10)


png("figures/GESMD_IWS_clustering/subgroup_prog_inter/RUNX1_aml_panel.png", res = 300, height = 2000, width = 2600)
plot_grid(
    plot_grid(ezh2_runx1_aml, stag2_runx1_aml, nrow = 1, labels = c("A", "B")),
     full_runx1_aml, ncol = 1, labels = c("", "C"))
dev.off()

## -7 and BM BLAST
mono7_blast_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "-7") %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), 
  variable = "BLAST", title = "BM BLAST in -7 subgroup", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts, pval_height = 0.10)

full_blast_aml <- makeSurvPlotAML(IWS_amlt %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), variable = "BLAST", 
  title = "BM BLAST in full cohort", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts, pval_height = 0.10)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/mono7_BMBLAST_aml_panel.png", res = 300, heigh = 1300, width = 2600)
plot_grid(mono7_blast_aml, full_blast_aml, ncol = 2)
dev.off()



























# plotPair <- function(group, gene){
#     df <- IWS_mds %>%
#         filter(sub_group == group) %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
    
#     p <- createSurvFit(df, gene)
#     df2 <- IWS_full %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
#     p2 <- createSurvFit(df2, gene)

#     df3 <- gesmd_dataset %>%
#         filter(sub_group == group) %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))   
#     p3 <- createSurvFit(df3, gene)

#     df4 <- gesmd_full %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
#     p4 <- createSurvFit(df4, gene)


#     plot_grid(
#         plot_grid(
#             makePanelPlot(p, paste(gene, "in", group, "- IWS")),
#             makePanelPlot(p2, paste(gene, "in all IWS")),
#             ncol = 2, rel_widths = c(0.8, 1)),
#         plot_grid(
#             makePanelPlot(p3, paste(gene, "in", group, "- GESMD")),
#             makePanelPlot(p4, paste(gene, "in all GESMD")),
#             ncol = 2, rel_widths = c(0.8, 1)),
#     ncol = 1)
# }


# plotPair <- function(group, gene, df_mds, df_full, title1, title2){
#     df <- df_mds %>%
#         filter(sub_group == group) %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
    
#     p <- createSurvFit(df, gene)
#     df2 <- df_full %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
#     p2 <- createSurvFit(df2, gene)

#     if (gene == "TET2other"){
#       gene = "TET2 mono-allelic"
#     }
#     if (gene == "plus8"){
#       gene = "8+"
#     }
#         plot_grid(
#             makePanelPlot(p, paste(gene, "in", group, title1)),
#             makePanelPlot(p2, paste(gene, title2)),
#             ncol = 2, rel_widths = c(0.8, 1))
# }



# computeModels <- function(var, group, outcome, dataset){

#   if(outcome == "OS"){
#     out = "OS_YEARS,OS_STATUS"
#   } else if (outcome == "AMLt"){
#     out = "AMLt_YEARS,AMLt_STATUS"
#   }

#   if (dataset == "IWS"){
#     df_mds = IWS_mds
#     df_full = IWS_full
#   } else if (dataset == "GESMD"){
#     df_mds = gesmd_dataset
#     df_full = gesmd_full
#   } else if (dataset == "joint"){
#     df_mds = joint_mds
#     df_full = joint_full
#   }
  
#   df_filt <- filter(df_mds, sub_group == group)

#   if (group == "7-"){
#     group = "del7"
#   }
#   if (group == "TET2 bi-allelic"){
#    group = "TET2bi"
# }
#   formula_base <- paste("Surv(", out, ") ~", var, "+ AGE + SEX", sep = "")
#   formula_interaction <- paste("Surv(", out, ") ~", var, "*", group, "+ AGE + SEX", sep = "")

#   if (dataset == "joint"){
#     formula_base <- paste(formula_base, "+ dataset", sep = "")
#     formula_interaction <- paste(formula_interaction, "+ dataset", sep = "")
#   }

#   model_subgroup <- coxph(formula(formula_base), df_filt)
#   model_full <- coxph(formula(formula_base), df_full)
#   model_interaction <- coxph(formula(formula_interaction), df_full)
#   return(list(subgroup = model_subgroup, full = model_full, interaction = model_interaction))
# }

# dataset <- c("IWS", "GESMD", "joint")
# names(dataset) <- dataset 

# ## EZH2 and TET2other
# ezh2_TET2other <- plotPair("EZH2", "TET2other", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_TET2other_joint.png", res = 300, heigh = 1500, width = 3000)
# ezh2_TET2other
# dev.off()

# ezh2_tet2other_models <- lapply(dataset, computeModels, var = "TET2other", group = "EZH2", outcome = "OS")

# # png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_TET2other.png", res = 300, heigh = 3000, width = 2800)
# # plotPair("EZH2", "TET2other", IWS_mds, IWS_full, "IWS", "IWS")
# # dev.off()


# ## STAG2 and RUNX1
# stag2_runx1_os <- plotPair("STAG2", "RUNX1", joint_mds, joint_full, "joint cohort", "MDS samples")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_RUNX1_joint.png", res = 300, heigh = 1400, width = 3000)
# stag2_runx1_os
# dev.off()

# stag2_runx1_models <- lapply(dataset, computeModels, var = "RUNX1", group = "STAG2", outcome = "OS")
# ## Significativo en IWS pero no en GESMD

# ## EZH2 and SETBP1
# ezh2_setbp1_os <- plotPair("EZH2", "SETBP1", joint_mds, joint_full, "joint cohort", "MDS samples")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_SETBP1_joint.png", res = 300, heigh = 1400, width = 2800)
# ezh2_setbp1_os
# dev.off()

# ezh2_setbp1_models <- lapply(dataset, computeModels, var = "SETBP1", group = "EZH2", outcome = "OS")

# ## STAG2 and U2AF1
# stag2_u2af1_os <- plotPair("STAG2", "U2AF1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_U2AF1_joint.png", res = 300, heigh = 1400, width = 2800)
# stag2_u2af1_os
# dev.off()

# stag2_u2af1_models <- lapply(dataset, computeModels, var = "U2AF1", group = "STAG2", outcome = "OS")

# ## EZH2 and U2AF1
# ezh2_u2af1_os <- plotPair("EZH2", "U2AF1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_U2AF1_joint.png", res = 300, heigh = 1400, width = 2800)
# ezh2_u2af1_os
# dev.off()

# ezh2_u2af1_models <- lapply(dataset, computeModels, var = "U2AF1", group = "EZH2", outcome = "OS")

# ## STAG2 and 8+
# stag2_plus8_os <- plotPair("STAG2", "plus8", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_plus8_joint.png", res = 300, heigh = 1400, width = 2800)
# stag2_plus8_os
# dev.off()

# stag2_plus8_models <- lapply(dataset, computeModels, var = "plus8", group = "STAG2", outcome = "OS")

# ## La interacción es significativa

# ## 7- and SETBP1
# del7_SETBP1_os <- plotPair("7-", "SETBP1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_SETBP1_joint.png", res = 300, heigh = 1400, width = 2800)
# del7_SETBP1_os
# dev.off()

# del7_setbp1_models <- lapply(dataset, computeModels, var = "SETBP1", group = "7-", outcome = "OS")
# ## Interacción significativa  

# ## U2AF1 and 7-
# del7_u2af1_os <- plotPair("7-", "U2AF1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_U2AF1_joint.png", res = 300, heigh = 1400, width = 2800)
# del7_u2af1_os
# dev.off()

# del7_u2af1_models <- lapply(dataset, computeModels, var = "U2AF1", group = "7-", outcome = "OS")

# ## TET2bi and ASXL1
# tet2bi_asxl1_os <- plotPair("TET2 bi-allelic", "ASXL1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/TET2bi_ASXL1_joint.png", res = 300, heigh = 1400, width = 2800)
# tet2bi_asxl1_os
# dev.off()

# tet2bi_asxl1_models <- lapply(dataset, computeModels, var = "ASXL1", group = "TET2 bi-allelic", outcome = "OS")

# ## Interacción no significativa

# ## 7- and ASXL1
# del7_asxl1_os <- plotPair("7-", "ASXL1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_ASXL1_joint.png", res = 300, heigh = 1400, width = 2800)
# del7_asxl1_os
# dev.off()

# del7_asxl1_models <- lapply(dataset, computeModels, var = "ASXL1", group = "7-", outcome = "OS")


# ## La interacción no es significativa

# ## EZH2 and RUNX1
# ezh2_runx1_os <- plotPair("EZH2", "RUNX1", joint_mds, joint_full, "joint cohort", "MDS samples")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_RUNX1_joint.png", res = 300, heigh = 1400, width = 2800)
# ezh2_runx1_os
# dev.off()

# ezh2_runx1_models <- lapply(dataset, computeModels, var = "RUNX1", group = "EZH2", outcome = "OS")

# ## CBL and EZH2
# ehz2_cbl_models <- lapply(dataset, computeModels, var = "CBL", group = "EZH2", outcome = "OS")

# ## TET2bi and RUNX1
# tet2bi_runx1_models <- lapply(dataset, computeModels, var = "RUNX1", group = "TET2 bi-allelic", outcome = "OS")

# ## TET2bi and U2AF1
# tet2bi_u2af1_models <- lapply(dataset, computeModels, var = "U2AF1", group = "TET2 bi-allelic", outcome = "OS")


# ## Clinical variables
# joint_full_plt <- mutate(joint_full,
#     PLT_cont = ifelse(dataset == "GESMD", exp(PLT), PLT),
#     PLT_cont = pmin(PLT_cont, 250, na.rm = TRUE)) %>%
#     filter(!is.na(PLT_cont))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*del7 + AGE + SEX, IWS_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*del7 + AGE + SEX + dataset, joint_full))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*STAG2 + AGE + SEX, IWS_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*STAG2 + AGE + SEX + dataset, joint_full))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*EZH2 + AGE + SEX, IWS_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*EZH2 + AGE + SEX + dataset, joint_full))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*TET2bi + AGE + SEX, IWS_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*TET2bi + AGE + SEX + dataset, joint_full))


# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ HB*del7 + AGE + SEX + dataset, joint_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ HB*STAG2 + AGE + SEX + dataset, joint_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ HB*EZH2 + AGE + SEX + dataset, joint_full))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ HB*TET2bi + AGE + SEX + dataset, joint_full))

# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ PLT_cont*del7 + AGE + SEX + dataset, joint_full_plt))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ PLT_cont*STAG2 + AGE + SEX + dataset, joint_full_plt))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ PLT_cont*EZH2 + AGE + SEX + dataset, joint_full_plt))
# summary(coxph(formula = Surv(OS_YEARS,OS_STATUS) ~ PLT_cont*TET2bi + AGE + SEX + dataset, joint_full_plt))

# ### STAG vs BM
# stag2_bm_models <- lapply(dataset, computeModels, var = "BM_BLAST", group = "STAG2", outcome = "OS")

# joint_mds_STAG2 <- filter(joint_mds, sub_group == "STAG2") %>%
#     mutate(BM = ifelse(BM_BLAST < 5, "<5%", 
#         ifelse(BM_BLAST < 10, "5-10%", "10%+")),
#         BM = factor(BM, levels = c("<5%", "5-10%", "10%+")))  %>%
#         filter(!is.na(BM) & !is.na(OS_YEARS) & !is.na(OS_STATUS))

# surv_STAG2 <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ BM, joint_mds_STAG2) %>%
#     ggsurvplot(data = joint_mds_STAG2, surv.median.line = "hv", 
#      risk.table = TRUE, break.time.by = 2, 
#       legend.labs  = levels(joint_mds_STAG2$BM))

# surv_STAG2_plot <- 
#     plot_grid(surv_STAG2$plot + 
#         ggtitle("BM BLASTS (%) in STAG2") +
#          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#          ylab("OS probability") +
#          xlab("Time (years)"), 
#         surv_STAG2$table, ncol = 1)

# png("figures/GESMD_IWS_clustering/subgroup_interaction/OS_STAG2_BM_joint.png", width = 4000, height = 2000, res = 300)
# surv_STAG2_plot
# dev.off()


# sel_inter_os <- mut_hr_comb_joint %>%
#     filter(Gene %in% c("TET2 mono-allelic", "RUNX1", "SETBP1", "U2AF1", "8+")) %>%
#   ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-Group") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Overall Survival in IWS")    

# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_sel_joint.png", res = 300, height = 1300, width = 2200)
# sel_inter_os
# dev.off()

# rest_inter_os <- mut_hr_comb_joint %>%
#     filter(!Gene %in% c("TET2 mono-allelic", "RUNX1", "SETBP1", "U2AF1", "8+")) %>%
#   ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-Group") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Overall Survival in IWS")    

# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_rest_joint.png", res = 300, heigh = 1300, width = 2200)
# rest_inter_os
# dev.off()


# os_signif_inter <- plot_grid(ezh2_TET2other, ezh2_setbp1_os,
#     stag2_runx1_os, stag2_u2af1_os, stag2_plus8_os, del7_SETBP1_os,
#     ncol = 2, labels = "AUTO")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/OS_surv_interactions_panel_joint.png", res = 300, height = 5000, width = 6000)
# os_signif_inter
# dev.off()

# ## AML transformation
# #####################################################################################
# mut_hr_aml_clust <- lapply(test_muts, function(gene){
#   lapply(groups2, function(cl){
#     df <- subset(IWS_mds, sub_group == cl)
#     df_est <- getEstimates(gene, df, outcome = "AMLt") %>%
#       mutate(sub_group = cl)
#   }) %>% Reduce(f = rbind)
# }) %>% Reduce(f = rbind)
# mut_hr_aml_full <- lapply(test_muts, function(gene){
#   df_est <- getEstimates(gene, IWS_full, outcome = "AMLt") %>%
#       mutate(sub_group = "IPSSM cohort")
#   }) %>% Reduce(f = rbind)

# mut_hr_aml_comb <- rbind(mut_hr_aml_clust, mut_hr_aml_full) %>%
#   mutate(sub_group = factor(sub_group, levels = c(groups2, "IPSSM cohort"))) %>%
#   filter(N >= 10) %>%
#   mutate(sub_group = droplevels(sub_group)) %>%
#   as_tibble() %>%
#   mutate(Gene = factor(case_when(
#     Gene == "plus8" ~ "8+",
#     Gene == "TET2other" ~ "TET2 mono-allelic",
#     TRUE ~ Gene
#   )))

# aml_int_plot <- mut_hr_aml_comb %>%
#     filter(HR > 0.1) %>%
#   ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-group") +
#   ggtitle("AML Transformation in IWS") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   scale_color_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5))    
# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_AMLt.png", res = 300, heigh = 1300, width = 2200)
#  aml_int_plot
# dev.off()


# mut_hr_aml_comb2 <- mut_hr_aml_comb %>%
#  group_by(Gene) %>% 
#  mutate(HR_L_Ref = HR_L[sub_group == "Low blasts"],
#         HR_H_Ref = HR_H[sub_group == "Low blasts"]) %>%
#  ungroup()
# top_low_aml <- filter(mut_hr_aml_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
#     arrange(HR_H - HR_L_Ref)
# top_high_aml <- filter(mut_hr_aml_comb2, HR < HR_L_Ref | HR > HR_H_Ref) %>%
#     arrange(HR_H_Ref - HR_L)
# filter(mut_hr_aml_comb2, HR_H < HR_L_Ref | HR_L > HR_H_Ref)


# createSurvFitAML <- function(df, gene){

#     mod <- summary(coxph(formula(paste("Surv(AMLt_YEARS,AMLt_STATUS) ~", gene, "+ AGE + SEX")), 
#                        df))
#     hr <- mod$coefficients[1, 2]
#     pval <- mod$coefficients[1, 5]
#     survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ gene, df) %>%
#         ggsurvplot(data = df, 
#                  risk.table = TRUE, break.time.by = 2, 
#                  fun = "event",
#                  pval = sprintf("HR = %.2f\nP = %.3f", hr, pval),
#                  pval.coord = c(max(df$AMLt_YEARS, na.rm = TRUE)*0.7, 0.2),
#                  legend.labs  = c("mut", "wt"), ylim = c(0, 1))  +
#         xlab("Time (Years)") 
# }
# makePanelPlot <- function(plot, title){
#     plot_grid(plot$plot + 
#                 theme(legend.position = "none", 
#                     plot.title = element_text(hjust = 0.5)) +
#                 ylab("AML transformation probability") +
#                 xlab("Time (years)") +
#                 ggtitle(title) +
#                 theme(plot.title = element_text(hjust = 0.5)), 
#                 plot$table, ncol = 1, rel_heights = c(1, 0.4))

# }

# plotPairAML <- function(group, gene){
#     df <- IWS_mds %>%
#         filter(sub_group == group) %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
    
#     p <- createSurvFitAML(df, gene)
#     df2 <- IWS_full %>%
#         mutate(gene = ifelse(get(gene) == 1, "mut", "WT")) %>%
#         filter(!is.na(gene))
#     p2 <- createSurvFitAML(df2, gene)
  
#   if (gene == "TET2other"){
#       gene = "TET2 mono-allelic"
#     }
#     if (gene == "plus8"){
#       gene = "8+"
#     }
#         plot_grid(
#             makePanelPlot(p, paste(gene, "in", group)),
#             makePanelPlot(p2, paste(gene, "in all IWS")),
#             ncol = 2, rel_widths = c(0.8, 1))
# }

# ## STAG2 and RUNX1
# stag2_runx1_aml <- plotPairAML("STAG2", "RUNX1")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_RUNX1_AMLt.png", res = 300, height = 1500, width = 2800)
# stag2_runx1_aml
# dev.off()

# stag2_runx1_aml_models <- computeModels(var = "RUNX1", group = "STAG2", outcome = "AMLt", dataset = "IWS")


# ## EZH2 and SETPBP1
# ezh2_setbp1_aml <- plotPairAML("EZH2", "SETBP1")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_SETBP1_AMLt.png", res = 300, height = 1500, width = 2800)
# ezh2_setbp1_aml
# dev.off()

# ezh2_setbp1_aml_models <- computeModels(var = "SETBP1", group = "EZH2", outcome = "AMLt", dataset = "IWS")

# ## STAG2 and U2AF1
# stag2_u2af1_aml <- plotPairAML("STAG2", "U2AF1")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_U2AF1_AMLt.png", res = 300, height = 1500, width = 2800)
# stag2_u2af1_aml
# dev.off()

# stag2_u2af1_aml_models <- computeModels(var = "U2AF1", group = "STAG2", outcome = "AMLt", dataset = "IWS")

# ## EZH2 and U2AF1
# ezh2_u2af1_aml <- plotPairAML("EZH2", "U2AF1")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_U2AF1_AMLt.png", res = 300, height = 1500, width = 2800)
# ezh2_u2af1_aml
# dev.off()

# ezh2_u2af1_aml_models <- computeModels(var = "U2AF1", group = "EZH2", outcome = "AMLt", dataset = "IWS")

# ## STAG and 8+
# stag2_plus8_aml <- plotPairAML("STAG2", "plus8")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/STAG2_dup8_AMLt.png", res = 300, height = 1500, width = 2800)
# stag2_plus8_aml
# dev.off()

# stag2_plus8_aml_models <- computeModels(var = "plus8", group = "STAG2", outcome = "AMLt", dataset = "IWS")

# ## SETBP1 and 7-
# del7_setbp1_aml <- plotPairAML("7-", "SETBP1")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_SETBP1_AMLt.png", res = 300, height = 1500, width = 2800)
# del7_setbp1_aml
# dev.off()

# del7_setbp1_aml_models <- computeModels(var = "SETBP1", group = "7-", outcome = "AMLt", dataset = "IWS")

# ## EZH2 and TET2other
# ezh2_TET2other_aml <- plotPairAML("EZH2", "TET2other")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_TET2other_AMLt.png", res = 300, height = 1500, width = 2800)
# ezh2_TET2other_aml
# dev.off()

# ezh2_tet2other_aml_models <- computeModels(var = "TET2other", group = "EZH2", outcome = "AMLt", dataset = "IWS")

# ## EZH2 and RUNX1 
# ezh2_runx1_aml <- plotPairAML("EZH2", "RUNX1")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/EZH2_RUNX1_AMLt.png", res = 300, height = 1500, width = 2800)
# ezh2_runx1_aml
# dev.off()

# ezh2_runx1_aml_models <- computeModels(var = "RUNX1", group = "EZH2", outcome = "AMLt", dataset = "IWS")

# ## ASXL1 and TET2bi
# asxl1_tet2bi_aml <- plotPairAML("TET2 bi-allelic", "ASXL1")
# png("figures/GESMD_IWS_clustering/subgroup_interaction/ASXL1_TET2bi_AMLt.png", res = 300, height = 1500, width = 2800)
# asxl1_tet2bi_aml
# dev.off()

# asxl1_tet2bi_aml_models <- computeModels(var = "ASXL1", group = "TET2 bi-allelic", outcome = "AMLt", dataset = "IWS")

# ## 7- and ASXL1
# del7_asxl1_aml <- plotPairAML("7-", "ASXL1")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/del7_ASXL1_AMLt.png", res = 300, height = 1500, width = 2800)
# del7_asxl1_aml
# dev.off()

# del7_asxl1_aml_models <- computeModels(var = "ASXL1", group = "7-", outcome = "AMLt", dataset = "IWS")

# ## STAG2 and ASXL1
# stag2_asxl1_aml_models <- computeModels(var = "ASXL1", group = "STAG2", outcome = "AMLt", dataset = "IWS")

# ## TET2bi and DNMT3A
# tet2bi_dnmt3a_aml_models <- computeModels(var = "DNMT3A", group = "TET2 bi-allelic", outcome = "AMLt", dataset = "IWS")

# ## 7- and DNMT3A (low N)
# del7_dnmt3a_aml_models <- computeModels(var = "DNMT3A", group = "7-", outcome = "AMLt", dataset = "IWS")

# ## Clinical variables

# ## BM
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ BM_BLAST*EZH2 + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ BM_BLAST*TET2bi + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ BM_BLAST*del7 + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ BM_BLAST*STAG2 + AGE + SEX, IWS_full)

# ## PLT 
# IWS_full_plt <- subset(joint_full_plt, dataset == "IWS")
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2*EZH2 + AGE + SEX, IWS_full_plt)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2*TET2bi + AGE + SEX, IWS_full_plt)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2*del7 + AGE + SEX, IWS_full_plt)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2*STAG2 + AGE + SEX, IWS_full_plt)

# ## HB
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ HB*EZH2 + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ HB*TET2bi + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ HB*del7 + AGE + SEX, IWS_full)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ HB*STAG2 + AGE + SEX, IWS_full)

# computeModels(var = "BM_BLAST", group = "STAG2", outcome = "AMLt", dataset = "IWS")
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2 + AGE + SEX, IWS_full_plt, subset = STAG2 == 1)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2 + AGE + SEX, IWS_full_plt)

# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2 + AGE + SEX, IWS_full_plt, subset = del7 == 1)
# coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT2 + AGE + SEX, IWS_full_plt)

# computeModels(var = "HB", group = "EZH2", outcome = "AMLt", dataset = "IWS")



# IWS_clin_models <- IWS_full_plt %>%
#     mutate(BM = ifelse(BM_BLAST < 5, "<5%", 
#         ifelse(BM_BLAST < 10, "5-10%", "10%+")),
#         BM = factor(BM, levels = c("<5%", "5-10%", "10%+")),
#         HB_cat = ifelse(HB < 10, "<10", 
#             ifelse(HB < 12, "10-12", "12+")),
#             HB_cat = factor(HB_cat, levels = c("<10", "10-12", "12+")),
#         PLT_cat = ifelse(PLT2 < 50, "<50", 
#             ifelse(PLT2 < 100, "50-100", 
#             ifelse(PLT2 > 150, "150+", "100-150"))),
#             PLT_cat = factor(PLT_cat, levels = c("<50", "50-100", "100-150", "150+"))) %>%
#       left_join(select(IWS_mds, ID, sub_group), by = "ID") 




# aml_STAG2_bm <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ BM, IWS_clin_models, subset = sub_group == "STAG2") %>%
#     ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#      risk.table = TRUE, break.time.by = 2, 
#       legend.labs  = levels(IWS_clin_models$BM))

# aml_STAG2_bm_plot <- 
#     plot_grid(aml_STAG2_bm$plot + 
#         ggtitle("BM BLASTS (%) in STAG2") +
#          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#          ylab("AML transformation") +
#          xlab("Time (years)"), 
#         aml_STAG2_bm$table, ncol = 1)

# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_STAG2_BM_joint.png", width = 4000, height = 2000, res = 300)
# aml_STAG2_bm_plot
# dev.off()


#  aml_STAG2_plt <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT_cat, IWS_clin_models, subset = sub_group == "STAG2") %>%
#    ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#     ylim= c(0, 1),
#      risk.table = TRUE, break.time.by = 2, 
#       legend.labs  = levels(IWS_clin_models$PLT_cat))

#  aml_STAG2_plt_full <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT_cat, IWS_clin_models) %>%
#    ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#     ylim = c(0, 1),
#      risk.table = TRUE, break.time.by = 2, 
#       legend.labs  = levels(IWS_clin_models$PLT_cat))


# aml_STAG2_plt_plot <- makePanelPlot(aml_STAG2_plt, "PLT in STAG2 samples")
# aml_STAG2_plt_full <- makePanelPlot(aml_STAG2_plt_full, "PLT in IWS cohort")

# aml_STAG2_plt_panel <- plot_grid(aml_STAG2_plt_plot, aml_STAG2_plt_full, ncol = 2, rel_widths = c(0.8, 1))


# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_STAG2_plt.png", width = 4000, height = 2000, res = 300)
# aml_STAG2_plt_panel
# dev.off()

# ## 7- and PLT
# aml_del7_plt <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ PLT_cat, IWS_clin_models, subset = sub_group == "7-") %>%
#    ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#      risk.table = TRUE, break.time.by = 2, 
#      ylim = c(0, 1),
#       legend.labs  = levels(IWS_clin_models$PLT_cat))

# aml_del7_plt_plot <- makePanelPlot(aml_del7_plt, "PLT in 7- samples")

# aml_plt_panel <- plot_grid(plot_grid(aml_STAG2_plt_plot, aml_del7_plt_plot, ncol = 2),
#  aml_del7_plt_plot_full, ncol = 1)

# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_plt.png", width = 4000, height = 4000, res = 300)
# aml_plt_panel
# dev.off()

# ## EZH2 and HB
# aml_ezh2_hb <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ HB_cat, IWS_clin_models, subset = sub_group == "EZH2") %>%
#    ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#      risk.table = TRUE, break.time.by = 2, 
#      ylim = c(0, 1),
#       legend.labs  = levels(IWS_clin_models$HB_cat))

# aml_ezh2_hb_full <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ HB_cat, IWS_clin_models) %>%
#    ggsurvplot(data = IWS_clin_models,
#     fun = "event",
#      risk.table = TRUE, break.time.by = 2, 
#      ylim = c(0, 1),
#       legend.labs  = levels(IWS_clin_models$HB_cat))

# aml_ezh2_hb_plot <- makePanelPlot(aml_ezh2_hb, "HB in EZH2 samples")
# aml_ezh2_hb_full_plot <- makePanelPlot(aml_ezh2_hb_full, "HB in IWS samples")

# aml_ezh2_hb_panel <- plot_grid(aml_ezh2_hb_plot, aml_ezh2_hb_full_plot, ncol = 2)
# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_EZH2_HB.png", width = 4000, height = 2000, res = 300)
# aml_ezh2_hb_panel
# dev.off()

# aml_clinical_panel <- plot_grid(aml_plt_panel,
#     aml_ezh2_hb_panel,
#     ncol = 1, labels = "AUTO", rel_heights = c(2, 1)
# )
# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_clinical_panel.png", width = 4000, height = 5000, res = 300)
# aml_clinical_panel
# dev.off()



# sel_inter_aml <- mut_hr_aml_comb %>%
#     filter(Gene %in% c("TET2 mono-allelic", "RUNX1", "8+")) %>%
#   ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-Group") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5)) +
#   ggtitle("AML transformation in IWS")    

# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_aml_sel.png", res = 300, height = 1300, width = 2200)
# sel_inter_aml
# dev.off()

# rest_inter_aml <- mut_hr_aml_comb %>%
#     filter(!Gene %in% c("TET2 mono-allelic", "RUNX1", "8+")) %>%
#     filter(HR > 0.1) %>%
#   ggplot(aes(x = sub_group, y = HR, fill = sub_group)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
#   theme_bw() +
#   scale_y_continuous(transform = "log2") +
#   xlab("Sub-Group") +
#   facet_wrap(~ Gene, scales = "free_x") +
#   scale_fill_manual(name = "", values = c(colors6, "#59758a")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#   plot.title = element_text(hjust = 0.5)) +
#   ggtitle("AML transformation in IWS")    

# png("figures/GESMD_IWS_clustering/subgroup_interaction/HR_mutations_aml_rest.png", res = 300, height = 1300, width = 2200)
# rest_inter_aml
# dev.off()


# aml_interaction_sel <- plot_grid(ezh2_TET2other_aml, ezh2_runx1_aml,
#     stag2_runx1_aml,  stag2_plus8_aml,
#     ncol = 2, labels = "AUTO")

# png("figures/GESMD_IWS_clustering/subgroup_interaction/AMLt_interactions_panel_joint.png", res = 300, height = 3000, width = 6000)
# aml_interaction_sel
# dev.off()


# png("figures/GESMD_IWS_clustering/subgroup_interaction/interaction_panel.png", res = 300, height = 4000, width = 4000)

# plot_grid(
#     plot_grid(sel_inter_os + theme(legend.position = "none"),
#         plotPairIWS("EZH2", "TET2other") + theme(legend.position = "none"),
#         plotPairIWS("STAG2", "RUNX1") + theme(legend.position = "none"),
#         ncol = 1, 
#         rel_heights = c(1, 1, 1), 
#         labels = "AUTO"),
#     plot_grid(sel_inter_aml,
#         plotPairAML("EZH2", "TET2other") + theme(legend.position = "none"),
#         plotPairAML("STAG2", "RUNX1") + theme(legend.position = "none"),
#         ncol = 1, 
#         rel_heights = c(1, 1, 1), 
#     labels = c("D","E","F")),
# ncol = 2)
# dev.off()



# prognosis_interaction_panel <- plot_grid(
#   plot_grid(
#     os_comb,
#     plot_grid(hr_os_group_plot2, hr_os_group_plot, ncol = 1, labels = c("B", "C")),
#     ncol = 2, labels = c("A", "")),
#   plot_grid(
#     aml_all,
#     plot_grid(hr_aml_group_plot2, hr_aml_group_plot, ncol = 1, labels = c("E", "F")),
#     ncol = 2, labels = c("D", "")
#   ),
#   plot_grid(
#     plot_grid(sel_inter_os + theme(legend.position = "none"), 
#       sel_inter_aml + theme(legend.position = "bottom"), 
#       ncol = 1, labels = c("G", "H"), rel_heights = c(2, 1.5)),
#     plot_grid(surv_STAG2_plot, aml_STAG2_bm_plot,
#       ncol = 1, labels = c("I", "J"))
#   ),
#   ncol = 1, rel_heights = c(1, 1, 1.4)
# ) 
# png("figures/GESMD_IWS_clustering/subgroup_interaction/prognosis_interaction_panel.png", res = 300, width = 6000, height = 6000)
# prognosis_interaction_panel
# dev.off()


# ## Modified score
# ipssm_process <- IPSSMprocess(gesmd_full)
# ipssm_res <- IPSSMmain(ipssm_process)
# ipssm_annot <- IPSSMannotate(ipssm_res)

# gesmd_full_IPSSM <- gesmd_full %>%
#     select(-IPSSM) %>%
#     left_join(ipssm_annot %>% 
#     mutate(IPSSM = IPSSMcat_mean,
#     IPSSM = gsub(" ", "-", IPSSM , fixed = TRUE),
#     IPSSM_SCORE = IPSSMscore_mean) %>%
#     select(ID, IPSSM, IPSSM_SCORE), by = "ID") 



# IWS_new_score <- IWS_full %>%
#  mutate(sub_group = classifySamples(.),
#   sub_group = ifelse(consensus %in% c("del5q", "mutated SF3B1", "Mutated TP53"), consensus, as.character(sub_group)),
#   IPSSM_score2 = ifelse(sub_group == "STAG2", 1,
#   ifelse(sub_group == "7-", IPSSM_SCORE + 0.5, IPSSM_SCORE)),
#   IPSSM2 = ifelse(sub_group == "STAG2", "High", as.character(IPSSM)),
#   IPSSM2 = factor(IPSSM2, levels = levels(IWS_full$IPSSM))) %>%
#   filter(!is.na(IPSSM_score2) & !is.na(IPSSM_SCORE) & !is.na(OS_YEARS) & !is.na(OS_STATUS))

# GESMD_new_score <- gesmd_full_IPSSM %>%
#  mutate(sub_group = classifySamples(.),
#   sub_group = ifelse(consensus %in% c("del5q", "mutated SF3B1", "Mutated TP53"), consensus, as.character(sub_group)),
#   IPSSM_score2 = ifelse(sub_group == "STAG2", 1,
#   ifelse(sub_group == "7-", IPSSM_SCORE + 0.5, IPSSM_SCORE)),
#   IPSSM2 = ifelse(sub_group == "STAG2", "High", as.character(IPSSM)),
#   IPSSM2 = factor(IPSSM2, levels = levels(IWS_full$IPSSM))) %>%
#   filter(!is.na(IPSSM_score2) & !is.na(IPSSM_SCORE) & !is.na(OS_YEARS) & !is.na(OS_STATUS))


# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE, IWS_new_score))
# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_score2, IWS_new_score))

# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE, GESMD_new_score))
# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_score2, GESMD_new_score))

# concordance(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE, IWS_new_score))
# concordance(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_score2, IWS_new_score))




# joint_mds_score <- subset(joint_new_score, ID %in% joint_prognosis$ID)


# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM, IWS_new_score))
# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM2, IWS_new_score))

# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM, GESMD_new_score))
# concordance(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM2, GESMD_new_score))




# joint_mds_score <- joint_prognosis %>%
#  mutate(IPSSM_score2 = ifelse(sub_group %in% c("EZH2", "STAG2"), 0.72, IPSSM_SCORE)) %>%
#   filter(!is.na(IPSSM_score2) & !is.na(IPSSM_SCORE) & !is.na(OS_YEARS) & !is.na(OS_STATUS))

# summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + dataset, joint_mds_score))
# summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_score2 + dataset, joint_mds_score))




# # ## Old code for prognosis models

# # ori_blasts <- filter(joint_prognosis, sub_group %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
# #   mutate(sub_group = droplevels(sub_group))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*IPSSM_SCORE + AGE + SEX + dataset, data = ori_blasts)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = ori_blasts)


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + AGE + SEX + dataset, data = joint_prognosis)


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = !IPSSM %in% c("Very-High"))


# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")





# # png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_dataset.png", width = 4000, height = 2000, res = 300)
# # survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ dataset, joint_prognosis) %>%
# #     ggsurvplot(data = joint_prognosis, surv.median.line = "hv",
# #      risk.table = TRUE, break.time.by = 2, 
# #       legend.labs  = levels(joint_prognosis$dataset))
# # dev.off()

# # IWS_os <- mutate(IWS_mds, sub_group = relevel(sub_group, ref = "Low blasts"), IPSSM = relevel(IPSSM, ref = "Low"))

# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low", "Moderate-Low"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low", "Moderate-Low", "Moderate-High"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = IPSSM %in% c("Moderate-Low", "Moderate-High"))

# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts" & IPSSM %in% c("High", "Very-High")))

# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts"))


# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group %in%  c("EZH2", "STAG2")))

# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2"))
# # IWS_os %>%
# #   filter(sub_group == "STAG2" & !IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))  %>%
# #   mutate(IPSSM = droplevels(IPSSM)) %>%
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = .)

# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts" & IPSSM %in% c("High", "Very-High")))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group %in% c("Low blasts", "EZH2", "STAG2") & IPSSM %in% c("High", "Very-High")))

# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "Low blasts"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "EZH2"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "STAG2"))
# # summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "TET2 bi-allelic"))

# # low <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts")
# # stag <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2")
# # tet2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE +  AGE + SEX, data = IWS_os, subset = sub_group == "TET2 bi-allelic")
# # ezh2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2")
# # del7 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "7-")
# # low_temp <- termplot(low, term = 1, se = TRUE, plot = FALSE)
# # stag_temp <- termplot(stag, term = 1, se = TRUE, plot = FALSE)
# # tet2_temp <- termplot(tet2, term = 1, se = TRUE, plot = FALSE)
# # ezh2_temp <- termplot(ezh2, term = 1, se = TRUE, plot = FALSE)
# # del7_temp <- termplot(del7, term = 1, se = TRUE, plot = FALSE)

# # rbind(low_temp$IPSSM_SCORE %>% mutate(mod = "low"), 
# #               stag_temp$IPSSM_SCORE %>% mutate(mod = "STAG2"),
# #               tet2_temp$IPSSM_SCORE %>% mutate(mod = "TET2 bi-allelic"),
# #               ezh2_temp$IPSSM_SCORE %>% mutate(mod = "EZH2"),
# #               del7_temp$IPSSM_SCORE %>% mutate(mod = "7-")) %>%
# # ggplot(aes(x = x, y = y, col = mod)) + geom_line() + geom_point() +
# # ggtitle("OS")


# # low_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts")
# # stag_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2")
# # tet2_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "TET2 bi-allelic")
# # low_temp_aml <- termplot(low_aml, term = 1, se = TRUE, plot = FALSE)
# # stag_temp_aml <- termplot(stag_aml, term = 1, se = TRUE, plot = FALSE)
# # tet2_temp_aml <- termplot(tet2_aml, term = 1, se = TRUE, plot = FALSE)

# # rbind(low_temp_aml$IPSSM_SCORE %>% mutate(mod = "low"), 
# #               stag_temp_aml$IPSSM_SCORE %>% mutate(mod = "STAG2"),
# #               tet2_temp_aml$IPSSM_SCORE %>% mutate(mod = "TET2 bi-allelic")) %>%
# # ggplot(aes(x = x, y = y, col = mod)) + geom_line() + geom_point() +
# # ggtitle("AMLt")
# # dev.off()




# # ## Survival by IPSSM (only IWS)
# # IPSSM_groups <- levels(IWS_mds$IPSSM)
# # names(IPSSM_groups) <- IPSSM_groups


# # survs_ipssm <- lapply(IPSSM_groups, function(cat){
# #   df <- filter(IWS_mds, IPSSM == cat)
# #   n_group <- df %>% 
# #     group_by(sub_group) %>%
# #     summarize(n = n())
# #   sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
# #   df <- filter(df, sub_group %in% sel_clusts)
# #   p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
# #     ggsurvplot(data = df, surv.median.line = "hv",
# #                palette = colors6[which(groups %in% sel_clusts)],
# #                risk.table = TRUE, break.time.by = 2) +
# #     xlab("Time (Years)")
# #   p
# # })

# # gesmd_dataset$IPSSM_combined <- case_when( 
# #     gesmd_dataset$IPSSM %in% c("Very-Low", "Low") ~ "VeryLow-Low",
# #     gesmd_dataset$IPSSM %in% c("High", "Very-High") ~ "High-VeryHigh",
# #     TRUE ~ gesmd_dataset$IPSSM
# # )

# # IPSSM_groups2 <- unique(gesmd_dataset$IPSSM_combined)
# # names(IPSSM_groups2) <- IPSSM_groups2

# # survs_ipssm_gesmd <- lapply(IPSSM_groups2, function(cat){
# #   df <- filter(gesmd_dataset, IPSSM_combined == cat)
# #   n_group <- df %>% 
# #     group_by(sub_group) %>%
# #     summarize(n = n())
# #   sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
# #   df <- filter(df, sub_group %in% sel_clusts)
# #   p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
# #     ggsurvplot(data = df, surv.median.line = "hv",
# #                palette = colors6[which(groups %in% sel_clusts)],
# #                risk.table = TRUE, break.time.by = 2) +
# #     xlab("Time (Years)")
# #   p
# # })

# # lapply(IPSSM_groups, function(ipssm){

# #         plot_grid(
# #             plot_grid(survs_ipssm[[ipssm]]$plot + 
# #                 ggtitle(paste(ipssm, "IWS")) +
# #                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
# #                 ylab("OS probability") +
# #                 xlab("Time (years)"), 
# #                 survs_ipssm[[ipssm]]$table, ncol = 1),
# #             ncol = 1)

              
# #     ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", ipssm, ".png"),
# #         width = 1800, height = 2000, dpi = 300, units = "px")
# # })


# # ##  Moderate low
# # joint_mod_low <- joint_prognosis %>%
# #     filter(IPSSM == "Moderate-Low") %>%
# #     filter(!sub_group %in% c("7-", "STAG2", "EZH2")) %>%
# #     mutate(sub_group = droplevels(sub_group))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_mod_low)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_mod_low %>% filter(dataset == "IWS"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_mod_low)


# # ## Moderate High
# # joint_mod_high <- joint_prognosis %>%
# #     filter(IPSSM == "Moderate-High") %>%
# #     filter(!sub_group %in% c("7-")) %>%
# #     mutate(sub_group = droplevels(sub_group))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_mod_high)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_mod_high %>% filter(dataset == "IWS"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_mod_high)

# # joint_high <- joint_prognosis %>%
# #     filter(IPSSM == "High") %>%
# #     mutate(sub_group = droplevels(sub_group))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_high)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_high %>% filter(dataset == "IWS"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_high)


# # joint_moderate <- joint_prognosis %>%
# #     filter(IPSSM %in% c("Moderate-Low", "Moderate-High")) %>%
# #     filter(!sub_group %in% c("7-")) %>%
# #     mutate(sub_group = droplevels(sub_group))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM + AGE + SEX + dataset, data = joint_moderate)
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_moderate %>% filter(dataset == "IWS"))
# # coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset +  AGE + SEX, data = joint_moderate)


# #   df <- filter(IWS_mds, IPSSM %in% c("Moderate-Low", "Moderate-High"))
# #   n_group <- df %>% 
# #     group_by(sub_group) %>%
# #     summarize(n = n())
# #   sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
# #   df <- filter(df, sub_group %in% sel_clusts)
# #   p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
# #     ggsurvplot(data = df, surv.median.line = "hv",
# #                palette = colors6[which(groups %in% sel_clusts)],
# #                risk.table = TRUE, break.time.by = 2) +
# #     xlab("Time (Years)")
# #   png(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_Moderate_IPSSM_IWS.png"),
# #       width = 1800, height = 2000, res = 300)

# #   plot_grid(
# #       plot_grid(p$plot + 
# #           ggtitle(paste("Moderate IPSS-M IWS")) +
# #           theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
# #           ylab("OS probability") +
# #           xlab("Time (years)"), 
# #           p$table, ncol = 1),
# #       ncol = 1)
# # dev.off()


# # ## Survival by subgroup 
# # ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")
# # ipssm_cols2 <- ipssm_cols[c(1, 3, 4, 6)]
# # survs_subgroups <- lapply(groups, function(group){
# #   df <- filter(IWS_mds, sub_group == group)
# #   n_ipssm <- df %>% 
# #     group_by(IPSSM) %>%
# #     summarize(n = n())
# #   sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM)
# #   df <- filter(df, IPSSM %in% sel_ipssm)
# #   p <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ IPSSM, df) %>%
# #     ggsurvplot(data = df, surv.median.line = "hv",
# #                palette = ipssm_cols[which(IPSSM_groups %in% sel_ipssm)],
# #                risk.table = TRUE, break.time.by = 2) +
# #     xlab("Time (Years)")
# #   p
# # })

# # survs_subgroups_gesmd <- lapply(groups[!groups %in% c("7-", "PHF6")], function(group){
# #   df <- filter(gesmd_dataset, sub_group == group)
# #   n_ipssm <- df %>% 
# #     group_by(IPSSM_combined) %>%
# #     summarize(n = n())
# #   sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM_combined)
# #   df <- filter(df, IPSSM_combined %in% sel_ipssm)
# #   p <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ IPSSM_combined, df) %>%
# #     ggsurvplot(data = df, surv.median.line = "hv",
# #                palette = ipssm_cols2[which(IPSSM_groups2 %in% sel_ipssm)],
# #                risk.table = TRUE, break.time.by = 2) +
# #     xlab("Time (Years)")
# #   p
# # })

# # lapply(names(survs_subgroups), function(group){
    
# #     plot_grid(
# #         plot_grid(survs_subgroups[[group]]$plot + 
# #             ggtitle(group) +
# #              theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
# #              ylab("OS probability") +
# #              xlab("Time (years)"), 
# #             survs_subgroups[[group]]$table, ncol = 1),
# #         ncol = 1)
# #     ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", group, ".png"),
# #         width = 2000, height = 2000, dpi = 300, units = "px")
# # })





# # new_test <- IWS_mds %>%
# #   mutate(groups = ifelse(sub_group %in% c("EZH2", "7-"), "High new", as.character(sub_group)),
# #   groups = relevel(as.factor(groups), ref = "High new"))

# # coxph(Surv(OS_YEARS, OS_STATUS) ~ groups + AGE + SEX, data = new_test)
# # coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups + AGE + SEX, data = new_test)
# # coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups + IPSSM_SCORE + AGE + SEX, data = new_test)
# # coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups*IPSSM_SCORE + AGE + SEX, data = new_test)

# # coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os)
# # coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = IWS_os)



# # IWS_mds2 <- mutate(IWS_mds, sub_group = relevel(sub_group, ref = "Low blasts"), 
# #   MOLECULAR_GROUP = factor(MOLECULAR_GROUP),
# #   MOLECULAR_GROUP = relevel(MOLECULAR_GROUP, ref = "No-event"),
# #   MOLECULAR_GROUP2 = factor(case_when(
# #     MOLECULAR_GROUP %in% c("\"-7/SETBP1\"", "EZH2-ASXL1", "bi-TET2", "IDH-STAG2") ~ MOLECULAR_GROUP,
# #     TRUE ~ consensus
# #   )), 
# #   MOLECULAR_GROUP2 = relevel(MOLECULAR_GROUP2, ref = "Low blasts"))
# # main_sub <- coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = IWS_mds2)
# # int_sub <- coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_mds2) 

# # main_mol <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP + IPSSM_SCORE + AGE + SEX, data = IWS_mds2) 
# # int_mol <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP*IPSSM_SCORE + AGE + SEX, data = IWS_mds2) 

# # main_mol2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP2 + IPSSM_SCORE + AGE + SEX, data = IWS_mds2) 
# # int_mol2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP2*IPSSM_SCORE + AGE + SEX, data = IWS_mds2) 


# # gesmd_dataset2 <- mutate(gesmd_dataset, sub_group = relevel(sub_group, ref = "Low blasts"), 
# #  MOLECULAR_GROUP2 = factor(case_when(
# #     mol_manual %in% c("7-/SETBP1", "EZH2/ASXL1", "TET2/SRSF2", "IDH/STAG2") ~ mol_manual,
# #     TRUE ~ consensus
# #   )), 
# #   MOLECULAR_GROUP2 = relevel(MOLECULAR_GROUP2, ref = "Low blasts"))


# # gesmd_sub <- coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = gesmd_dataset2)
# # gesmd_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = gesmd_dataset2) 

# # gesmd_mol_sub <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP2 + IPSSM_SCORE + AGE + SEX, data = gesmd_dataset2) 
# # gesmd_mol_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ MOLECULAR_GROUP2*IPSSM_SCORE + AGE + SEX, data = gesmd_dataset2) 
