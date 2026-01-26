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

load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

## Modify reference levels
#IWS_mds$sub_group <- relevel(IWS_mds$sub_group, ref = "Low blasts")


colors6 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#999999", "grey40",  "black")

surv_IWS <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, IWS_mds) %>%
    ggsurvplot(data = IWS_mds, surv.median.line = "hv", palette = colors6,
     risk.table = TRUE, break.time.by = 2, 
      legend.labs  = levels(IWS_mds$sub_group))

surv_gesmd <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ sub_group, gesmd_dataset) %>%
    ggsurvplot(data = gesmd_dataset, surv.median.line = "hv", 
               palette = colors6, risk.table = TRUE, break.time.by = 2, 
      legend.labs  = levels(gesmd_dataset$sub_group)) 

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
    ncol = 2, rel_widths = c(1, 0.8))

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_all_subgroups.png", width = 4000, height = 2000, res = 300)
os_comb
dev.off()

joint_prognosis_plot <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), IPSSM, IPSSM_SCORE,  sub_group, AGE, SEX, dataset, BM_BLAST),
    gesmd_dataset %>% mutate(dataset = "GESMD") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), IPSSM, IPSSM_SCORE,  sub_group, AGE, SEX, dataset, BM_BLAST)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
    group2 = case_when(
        sub_group %in% c("Low blasts", "MDS-IB1", "MDS-IB2") ~ "Original",
        TRUE ~ sub_group
    ) %>% factor(levels = c("Original", c("EZH2", "TET2 bi-allelic", "7-", "STAG2"))),
     IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 


joint_prognosis <- joint_prognosis_plot %>% 
  mutate(sub_group = relevel(sub_group, ref = "Low blasts"), 
  IPSSM = factor(IPSSM, levels = c("Low", "Very-Low", "Moderate-Low", "Moderate-High", "High", "Very-High")))


raw <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group, data = joint_prognosis))
mod1 <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_prognosis))
main_os_sum <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) )
int_os_sum <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) )

makeModelTab <- function(cox_sum){

  coefs <- cox_sum$coefficients
  conf_int <- cox_sum$conf.int

  tab <- tibble(
    Variable = gsub("sub_group", "", rownames(coefs)),
    HR = sprintf("%.1f (%.1f - %.1f)", 
      round(coefs[, "exp(coef)"], 2), round(conf_int[, "lower .95"], 2), round(conf_int[, "upper .95"], 2)),
    p_value = signif(coefs[, "Pr(>|z|)"], 2)
  )
  tab
}

mods <- list(raw = raw, mod1 = mod1, main_os_sum = main_os_sum, int_os_sum = int_os_sum)
os_res <- lapply(names(mods), function(x) makeModelTab(mods[[x]]) %>% mutate(model = x)) %>% Reduce(f = rbind)

groups <- levels(joint_prognosis$sub_group)
names(groups) <- groups

os_res_filt <- os_res %>% 
  filter(str_detect(Variable, paste(groups, collapse = "|")))
write.table(os_res_filt, 
            file = "results/GESMD_IWS_clustering/OS_subgroups_coxph.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = TRUE,
            row.names = FALSE)

main_model_tabs <- lapply(sub_groups, function(group){

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
  tibble(Variable = "Moderate-High", Group = sub_groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
ggplot(aes(x = Variable, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors6) +
  geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "IPSSM",
    color = "Sub-group") +
  ggtitle("OS (Ref: IPSSM Moderate-High)")

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_os_group_plot
dev.off()


joint_prognosis_os_plot <- joint_prognosis %>%
    mutate(comb_group = paste(sub_group, IPSSM, sep = "_"),
    comb_group = factor(comb_group, levels = unique(comb_group)),
    comb_group = relevel(comb_group, ref = "Low blasts_Low"))

cat_count <- table(joint_prognosis_os_plot$comb_group)
sel_cats <- names(cat_count[cat_count >= 10])

joint_prognosis_os_plot <- joint_prognosis_os_plot %>%
  filter(comb_group %in% sel_cats)  %>%
  mutate(comb_group = droplevels(comb_group))
os_mod_subgroup_ipssm <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ comb_group + AGE + SEX + dataset, data = joint_prognosis_os_plot))

os_mod_subgroup_ipssm_tab <- tibble(
    Variable = gsub("comb_group", "", rownames(os_mod_subgroup_ipssm$coefficients)),
    Group = sapply(strsplit(Variable, "_", 2), `[`, 1),
    IPSSM = sapply(strsplit(Variable, "_", 2), `[`, 2),
    HR = round(os_mod_subgroup_ipssm$coefficients[, "exp(coef)"], 2),
    HR_Low = round(os_mod_subgroup_ipssm$conf.int[, "lower .95"], 2), 
    HR_High = round(os_mod_subgroup_ipssm$conf.int[, "upper .95"], 2),
    p_value = signif(os_mod_subgroup_ipssm$coefficients[, "Pr(>|z|)"], 2)
  ) %>%
  filter(IPSSM != "NA") %>%
  filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
  bind_rows(., tibble(
    Variable = "Low blasts_Low",
    Group = "Low blasts",
    IPSSM = "Low",
    HR = 1,
    HR_Low = 1,
    HR_High = 1,
    p_value = NA
  ))



hr_os_group_plot2 <- os_mod_subgroup_ipssm_tab %>%
mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
ggplot(aes(x = Group, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors6) +
  geom_errorbar(aes(x = Group, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ IPSSM, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sub-group",
  color = "Sub-group") +
  ggtitle("OS (Ref: Low blasts/IPSSM Low)")

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
hr_os_group_plot2
dev.off()

## OS by IPSSM (only IWS)
IPSSM_groups <- levels(IWS_mds$IPSSM)
names(IPSSM_groups) <- IPSSM_groups

groups <- levels(IWS_mds$sub_group)
os_joint_IPSSM <- lapply(IPSSM_groups, function(cat){
  df <- filter(joint_prognosis_plot, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors6[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
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

              
    ggsave(plot = os_plot, filename = paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", ipssm, "_joint.png"),
        width = 1800, height = 2000, dpi = 300, units = "px")
    os_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_IPSSM_joint_panel.png", width = 4000, height = 6000, res = 300)
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
               risk.table = TRUE, break.time.by = 2) +
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
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", group, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
    os_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_panel.png", width = 4000, height = 6000, res = 300)
plot_grid(plotlist = survs_subgroups_plots, ncol = 2, labels = "AUTO")
dev.off()

## AML transformation
amlt_IWS <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group, IWS_mds) %>%
    ggsurvplot(data = IWS_mds, surv.median.line = "hv", palette = colors6,
     risk.table = TRUE, break.time.by = 2, fun = "event",
      legend.labs  = levels(IWS_mds$sub_group))

aml_all <- plot_grid(amlt_IWS$plot + 
        ggtitle("IWS") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("AMLt probability") +
         xlab("Time (years)"), 
        amlt_IWS$table, ncol = 1)
png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_all_subgroups.png", width = 2000, height = 2000, res = 300)
aml_all
dev.off()

IWS_aml <- mutate(IWS_mds, sub_group = relevel(sub_group, ref = "Low blasts"), IPSSM = relevel(IPSSM, ref = "Low"))

aml_raw <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group, data = IWS_aml))
aml_adj <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group + AGE + SEX, data = IWS_aml))
aml_main <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = IWS_aml) )
aml_int <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX , data = IWS_aml) )


aml_mods <- list(raw = aml_raw, adj = aml_adj, main_aml = aml_main, int_aml = aml_int)
aml_res <- lapply(names(aml_mods), function(x) makeModelTab(aml_mods[[x]]) %>% mutate(model = x)) %>% Reduce(f = rbind)

aml_res_filt <- aml_res %>% 
  filter(str_detect(Variable, paste(groups, collapse = "|")))
write.table(aml_res_filt, 
            file = "results/GESMD_IWS_clustering/AMLt_subgroups_coxph.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = TRUE,
            row.names = FALSE)

aml_model_tabs <- lapply(sub_groups, function(group){

  IWS_sub <- IWS_aml %>%
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
  tibble(Variable = "Moderate-High", Group = sub_groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
      filter(HR_Low != 0) %>%
      filter(Group != "EZH2") %>%
ggplot(aes(x = Variable, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors6[-1]) +
  geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "IPSSM",
    color = "Sub-group") +
  ggtitle("AMLt (Reference: IPSSM Moderate-High)")

png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_aml_group_plot
dev.off()


IWS_aml_plot <- IWS_mds %>%
    mutate(comb_group = paste(sub_group, IPSSM, sep = "_"),
    comb_group = factor(comb_group, levels = unique(comb_group)),
    comb_group = relevel(comb_group, ref = "Low blasts_Low"))

cat_count_aml <- table(IWS_aml_plot$comb_group)
sel_cats_aml <- names(cat_count_aml[cat_count_aml >= 10])
sel_cats_aml <- sel_cats_aml[!sel_cats_aml %in% c("EZH2_Moderate-High", "Low blasts_NA", "TET2 bi-allelic_Very-Low")]
IWS_aml_plot <- IWS_aml_plot %>%
  filter(comb_group %in% sel_cats_aml)  %>%
  mutate(comb_group = droplevels(comb_group))
aml_mod_subgroup_ipssm <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ comb_group + AGE + SEX, data = IWS_aml_plot))

aml_mod_subgroup_ipssm_tab <- tibble(
    Variable = gsub("comb_group", "", rownames(aml_mod_subgroup_ipssm$coefficients)),
    Group = sapply(strsplit(Variable, "_", 2), `[`, 1),
    IPSSM = sapply(strsplit(Variable, "_", 2), `[`, 2),
    HR = round(aml_mod_subgroup_ipssm$coefficients[, "exp(coef)"], 2),
    HR_Low = round(aml_mod_subgroup_ipssm$conf.int[, "lower .95"], 2), 
    HR_High = round(aml_mod_subgroup_ipssm$conf.int[, "upper .95"], 2),
    p_value = signif(aml_mod_subgroup_ipssm$coefficients[, "Pr(>|z|)"], 2)
  ) %>%
  filter(IPSSM != "NA") %>%
  filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
  bind_rows(., tibble(
    Variable = "Low blasts_Low",
    Group = "Low blasts",
    IPSSM = "Low",
    HR = 1,
    HR_Low = 1,
    HR_High = 1,
    p_value = NA
  ))


hr_aml_group_plot2 <- aml_mod_subgroup_ipssm_tab %>%
mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
      filter(!IPSSM == "Very-Low") %>%
ggplot(aes(x = Group, y = HR, color = Group)) +
  geom_point() +
  scale_color_manual(values = colors6) +
  geom_errorbar(aes(x = Group, ymin = HR_Low, ymax = HR_High)) +
  theme_bw() +
  facet_grid(. ~ IPSSM, scale = "free_x", space = "free_x") +
  scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sub-group",
  color = "Sub-group") +
  ggtitle("AMLt (Ref: Low blasts/IPSSM Low)")

png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
hr_aml_group_plot2
dev.off()

panel_prognosis <- 
plot_grid(
  plot_grid(
    os_comb,
    plot_grid(hr_os_group_plot2, hr_os_group_plot, ncol = 1, labels = c("B", "C")),
    ncol = 2, labels = c("A", "")),
  plot_grid(
    aml_all,
    plot_grid(hr_aml_group_plot2, hr_aml_group_plot, ncol = 1, labels = c("E", "F")),
    ncol = 2, labels = c("D", "")
  ),
  ncol = 1, labels = "")

png("figures/GESMD_IWS_clustering/subgroup_prognosis/prognosis_panel.png", width = 4500, height = 4000, res = 300)
panel_prognosis
dev.off()


## AMLt by IPSSM (only IWS)
amlt_ipssm <- lapply(IPSSM_groups, function(cat){
  df <- filter(IWS_mds, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
               palette = colors6[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
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
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_", ipssm, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
    aml_plot
})

png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_IPSSM_panel.png", width = 4000, height = 6000, res = 300)
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








ori_blasts <- filter(joint_prognosis, sub_group %in% c("Low blasts", "MDS-IB1", "MDS-IB2")) %>%
  mutate(sub_group = droplevels(sub_group))
coxph(Surv(OS_YEARS,OS_STATUS) ~ BM_BLAST*IPSSM_SCORE + AGE + SEX + dataset, data = ori_blasts)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = ori_blasts)


coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + AGE + SEX + dataset, data = joint_prognosis)


coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2 + IPSSM + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")


coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = !IPSSM %in% c("Very-High"))


coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis)
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "IWS")
coxph(Surv(OS_YEARS,OS_STATUS) ~ group2*IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = dataset == "GESMD")





png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_dataset.png", width = 4000, height = 2000, res = 300)
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ dataset, joint_prognosis) %>%
    ggsurvplot(data = joint_prognosis, surv.median.line = "hv",
     risk.table = TRUE, break.time.by = 2, 
      legend.labs  = levels(joint_prognosis$dataset))
dev.off()

IWS_os <- mutate(IWS_mds, sub_group = relevel(sub_group, ref = "Low blasts"), IPSSM = relevel(IPSSM, ref = "Low"))

coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low", "Moderate-Low"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = !IPSSM %in% c("Low", "Very-Low", "Moderate-Low", "Moderate-High"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = IPSSM %in% c("Moderate-Low", "Moderate-High"))

summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts" & IPSSM %in% c("High", "Very-High")))

summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts"))


summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group %in%  c("EZH2", "STAG2")))

summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2"))
IWS_os %>%
  filter(sub_group == "STAG2" & !IPSSM %in% c("Very-Low", "Low", "Moderate-Low"))  %>%
  mutate(IPSSM = droplevels(IPSSM)) %>%
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = .)

summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts" & IPSSM %in% c("High", "Very-High")))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group %in% c("Low blasts", "EZH2", "STAG2") & IPSSM %in% c("High", "Very-High")))

summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "Low blasts"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "EZH2"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "STAG2"))
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "TET2 bi-allelic"))

low <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts")
stag <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2")
tet2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE +  AGE + SEX, data = IWS_os, subset = sub_group == "TET2 bi-allelic")
ezh2 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "EZH2")
del7 <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "7-")
low_temp <- termplot(low, term = 1, se = TRUE, plot = FALSE)
stag_temp <- termplot(stag, term = 1, se = TRUE, plot = FALSE)
tet2_temp <- termplot(tet2, term = 1, se = TRUE, plot = FALSE)
ezh2_temp <- termplot(ezh2, term = 1, se = TRUE, plot = FALSE)
del7_temp <- termplot(del7, term = 1, se = TRUE, plot = FALSE)

rbind(low_temp$IPSSM_SCORE %>% mutate(mod = "low"), 
              stag_temp$IPSSM_SCORE %>% mutate(mod = "STAG2"),
              tet2_temp$IPSSM_SCORE %>% mutate(mod = "TET2 bi-allelic"),
              ezh2_temp$IPSSM_SCORE %>% mutate(mod = "EZH2"),
              del7_temp$IPSSM_SCORE %>% mutate(mod = "7-")) %>%
ggplot(aes(x = x, y = y, col = mod)) + geom_line() + geom_point() +
ggtitle("OS")


low_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "Low blasts")
stag_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "STAG2")
tet2_aml <- coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = IWS_os, subset = sub_group == "TET2 bi-allelic")
low_temp_aml <- termplot(low_aml, term = 1, se = TRUE, plot = FALSE)
stag_temp_aml <- termplot(stag_aml, term = 1, se = TRUE, plot = FALSE)
tet2_temp_aml <- termplot(tet2_aml, term = 1, se = TRUE, plot = FALSE)

rbind(low_temp_aml$IPSSM_SCORE %>% mutate(mod = "low"), 
              stag_temp_aml$IPSSM_SCORE %>% mutate(mod = "STAG2"),
              tet2_temp_aml$IPSSM_SCORE %>% mutate(mod = "TET2 bi-allelic")) %>%
ggplot(aes(x = x, y = y, col = mod)) + geom_line() + geom_point() +
ggtitle("AMLt")
dev.off()




## Survival by IPSSM (only IWS)
IPSSM_groups <- levels(IWS_mds$IPSSM)
names(IPSSM_groups) <- IPSSM_groups


survs_ipssm <- lapply(IPSSM_groups, function(cat){
  df <- filter(IWS_mds, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors6[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})

gesmd_dataset$IPSSM_combined <- case_when( 
    gesmd_dataset$IPSSM %in% c("Very-Low", "Low") ~ "VeryLow-Low",
    gesmd_dataset$IPSSM %in% c("High", "Very-High") ~ "High-VeryHigh",
    TRUE ~ gesmd_dataset$IPSSM
)

IPSSM_groups2 <- unique(gesmd_dataset$IPSSM_combined)
names(IPSSM_groups2) <- IPSSM_groups2

survs_ipssm_gesmd <- lapply(IPSSM_groups2, function(cat){
  df <- filter(gesmd_dataset, IPSSM_combined == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors6[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})

lapply(IPSSM_groups, function(ipssm){

        plot_grid(
            plot_grid(survs_ipssm[[ipssm]]$plot + 
                ggtitle(paste(ipssm, "IWS")) +
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                ylab("OS probability") +
                xlab("Time (years)"), 
                survs_ipssm[[ipssm]]$table, ncol = 1),
            ncol = 1)

              
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", ipssm, ".png"),
        width = 1800, height = 2000, dpi = 300, units = "px")
})


##  Moderate low
joint_mod_low <- joint_prognosis %>%
    filter(IPSSM == "Moderate-Low") %>%
    filter(!sub_group %in% c("7-", "STAG2", "EZH2")) %>%
    mutate(sub_group = droplevels(sub_group))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_mod_low)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_mod_low %>% filter(dataset == "IWS"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_mod_low)


## Moderate High
joint_mod_high <- joint_prognosis %>%
    filter(IPSSM == "Moderate-High") %>%
    filter(!sub_group %in% c("7-")) %>%
    mutate(sub_group = droplevels(sub_group))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_mod_high)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_mod_high %>% filter(dataset == "IWS"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_mod_high)

joint_high <- joint_prognosis %>%
    filter(IPSSM == "High") %>%
    mutate(sub_group = droplevels(sub_group))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX + dataset, data = joint_high)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_high %>% filter(dataset == "IWS"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset + AGE + SEX, data = joint_high)


joint_moderate <- joint_prognosis %>%
    filter(IPSSM %in% c("Moderate-Low", "Moderate-High")) %>%
    filter(!sub_group %in% c("7-")) %>%
    mutate(sub_group = droplevels(sub_group))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + IPSSM + AGE + SEX + dataset, data = joint_moderate)
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group + AGE + SEX, data = joint_moderate %>% filter(dataset == "IWS"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*dataset +  AGE + SEX, data = joint_moderate)


  df <- filter(IWS_mds, IPSSM %in% c("Moderate-Low", "Moderate-High"))
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors6[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  png(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_Moderate_IPSSM_IWS.png"),
      width = 1800, height = 2000, res = 300)

  plot_grid(
      plot_grid(p$plot + 
          ggtitle(paste("Moderate IPSS-M IWS")) +
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
          ylab("OS probability") +
          xlab("Time (years)"), 
          p$table, ncol = 1),
      ncol = 1)
dev.off()


## Survival by subgroup 
ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")
ipssm_cols2 <- ipssm_cols[c(1, 3, 4, 6)]
survs_subgroups <- lapply(groups, function(group){
  df <- filter(IWS_mds, sub_group == group)
  n_ipssm <- df %>% 
    group_by(IPSSM) %>%
    summarize(n = n())
  sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM)
  df <- filter(df, IPSSM %in% sel_ipssm)
  p <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ IPSSM, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = ipssm_cols[which(IPSSM_groups %in% sel_ipssm)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})

survs_subgroups_gesmd <- lapply(groups[!groups %in% c("7-", "PHF6")], function(group){
  df <- filter(gesmd_dataset, sub_group == group)
  n_ipssm <- df %>% 
    group_by(IPSSM_combined) %>%
    summarize(n = n())
  sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM_combined)
  df <- filter(df, IPSSM_combined %in% sel_ipssm)
  p <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ IPSSM_combined, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = ipssm_cols2[which(IPSSM_groups2 %in% sel_ipssm)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})

lapply(names(survs_subgroups), function(group){
    
    plot_grid(
        plot_grid(survs_subgroups[[group]]$plot + 
            ggtitle(group) +
             theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
             ylab("OS probability") +
             xlab("Time (years)"), 
            survs_subgroups[[group]]$table, ncol = 1),
        ncol = 1)
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", group, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
})





new_test <- IWS_mds %>%
  mutate(groups = ifelse(sub_group %in% c("EZH2", "7-"), "High new", as.character(sub_group)),
  groups = relevel(as.factor(groups), ref = "High new"))

coxph(Surv(OS_YEARS, OS_STATUS) ~ groups + AGE + SEX, data = new_test)
coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups + AGE + SEX, data = new_test)
coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups + IPSSM_SCORE + AGE + SEX, data = new_test)
coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ groups*IPSSM_SCORE + AGE + SEX, data = new_test)

coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX, data = IWS_os)
coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group + IPSSM_SCORE + AGE + SEX, data = IWS_os)

