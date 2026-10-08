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
library(ggh4x)


load("results/GESMD_IWS_clustering/gesmd_IWS_mds.Rdata")
load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/hershberger/hershberger_mds.Rdata")
load("results/hershberger/hershberger_full.Rdata")


## Remove 3 classical molecular groups
clusters <- c("TET2-bi",  "-7", "EZH2", "STAG2", "MDS-LB", "MDS-IB1", "MDS-IB2")

IWS_mds_f <- IWS_mds %>% 
  filter(sub_group %in% clusters) %>%
  mutate(sub_group = fct_relevel(droplevels(sub_group), clusters))
gesmd_dataset_f <- gesmd_dataset %>% filter(sub_group %in% clusters) %>%
  mutate(sub_group = fct_relevel(droplevels(sub_group), clusters))
hersh_mds_f <- hersh_mds %>% filter(sub_group %in% clusters) %>%
  mutate(sub_group = factor(sub_group, levels = clusters))




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
# colors_all <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
#     "#D55E00", "#999999", "grey40",  "black")
colors <-  c("#56B4E9", "#009E73", "#E69F00",  "#CC79A7", 
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
    select(ends_with("STATUS"), ends_with("YEARS"), STAG2, IPSSM, IPSSM_SCORE,  sub_group, AGE, SEX, dataset, BM_BLAST),
    gesmd_dataset_f %>% mutate(dataset = "GESMD") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), STAG2, IPSSM, IPSSM_SCORE,  sub_group, AGE, SEX, dataset, BM_BLAST)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
     IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 


joint_prognosis <- joint_prognosis_plot %>% 
  mutate(IPSSM = factor(IPSSM, levels = c("Low", "Very-Low", "Moderate-Low", "Moderate-High", "High", "Very-High")))

groups <- levels(joint_prognosis$sub_group)
names(groups) <- groups

# main_model_tabs <- lapply(groups, function(group){

#   joint_prognosis_sub <- joint_prognosis %>%
#     filter(sub_group == group) %>%
#     mutate(sub_group = droplevels(sub_group),
#     IPSSM = relevel(IPSSM, ref = "Moderate-High"))
#     mod <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX + dataset, data = joint_prognosis_sub))
  
#   cat_sum <- table(joint_prognosis_sub$IPSSM)
#   sel_cats <- cat_sum[cat_sum >= 10] %>% names()
#   coefs <- mod$coefficients
#   conf_int <- mod$conf.int
#   tab <- tibble(
#     Variable = gsub("IPSSM", "", rownames(coefs)),
#     Group = group,
#     HR = round(coefs[, "exp(coef)"], 2),
#     HR_Low = round(conf_int[, "lower .95"], 2), 
#     HR_High = round(conf_int[, "upper .95"], 2),
#     p_value = signif(coefs[, "Pr(>|z|)"], 2)
#   ) %>%
#   filter(!Variable %in% c("AGE", "SEXM", "datasetGESMD")) %>%
#   filter(Variable %in% sel_cats)
#   tab
# }) %>% Reduce(f = rbind)


# hr_os_group_plot <- bind_rows(main_model_tabs,
#   tibble(Variable = "Moderate-High", Group = groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
# )  %>%
# mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
#       Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
# ggplot(aes(x = Variable, y = HR, color = Group)) +
#   geom_point() +
#   scale_color_manual(values = colors) +
#   geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
#   theme_bw() +
#   facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
#   scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#   plot.title = element_text(hjust = 0.5)) +
#   labs(x = "IPSSM",
#     color = "Sub-group") +
#   ggtitle("OS (Ref: IPSSM Moderate-High)")

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
# hr_os_group_plot
# dev.off()


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


median_OS_plot <- ggplot(ipssm_raw_surv_tib, aes(y = IPSSM, x = median, color = sub_group)) +
  geom_point(size = 2) +
  geom_segment(aes(x = `X0.95LCL`, xend = UCL),
  arrow = arrow(length = unit(ifelse(is.na(ipssm_raw_surv_tib$`X0.95UCL`), 0.3, 0), "cm"))) +
  scale_color_manual(values = colors) +
  labs(title = "Median OS by sub-group and IPSSM",
       y = "", color = "",
       x = "Median OS (years)") +
  facet_grid2(sub_group ~ ., scales = "free_y", 
  switch = "y", 
  space = "free_y",
    strip = strip_themed(
                        background_y = elem_list_rect(fill = colors),
                         text_y = element_text(angle = 0, color = "white", face = "bold")
                )) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
        scale_y_discrete(name = "", position = "right")

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
               palette = colors[which(clusters %in% sel_clusts)],
               legend.labs = clusters[which(clusters %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2, xlim = c(0, 10),
               xlab = "Time (Years)") 
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
              
    ggsave(plot = os_plot, filename = paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/split/OS_subgroups_", ipssm, "_joint.png"),
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
               risk.table = TRUE, break.time.by = 2, xlim = c(0, 10), xlab = "Time (Years)")
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
    # ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/OS_subgroups_", group, ".png"),
    #     width = 2000, height = 2000, dpi = 300, units = "px")
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

# aml_model_tabs <- lapply(groups, function(group){

#   IWS_sub <- IWS_mds_f %>%
#     filter(sub_group == group) %>%
#     mutate(sub_group = droplevels(sub_group),
#     IPSSM = relevel(IPSSM, ref = "Moderate-High"))
#     mod <- summary(coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ IPSSM + AGE + SEX, data = IWS_sub))
  
#   cat_sum <- table(IWS_sub$IPSSM)
#   sel_cats <- cat_sum[cat_sum >= 10] %>% names()
#   coefs <- mod$coefficients
#   conf_int <- mod$conf.int
#   tab <- tibble(
#     Variable = gsub("IPSSM", "", rownames(coefs)),
#     Group = group,
#     HR = round(coefs[, "exp(coef)"], 2),
#     HR_Low = round(conf_int[, "lower .95"], 2), 
#     HR_High = round(conf_int[, "upper .95"], 2),
#     p_value = signif(coefs[, "Pr(>|z|)"], 2)
#   ) %>%
#   filter(!Variable %in% c("AGE", "SEXM")) %>%
#   filter(Variable %in% sel_cats)
#   tab
# }) %>% Reduce(f = rbind)


# hr_aml_group_plot <- bind_rows(aml_model_tabs,
#   tibble(Variable = "Moderate-High", Group = groups, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
# )  %>%
# mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
#       Group = factor(Group, levels = levels(joint_prognosis_plot$sub_group))) %>%
#       filter(HR_Low != 0) %>%
#       filter(Group != "EZH2") %>%
# ggplot(aes(x = Variable, y = HR, color = Group)) +
#   geom_point() +
#   scale_color_manual(values = colors[-1]) +
#   geom_errorbar(aes(x = Variable, ymin = HR_Low, ymax = HR_High)) +
#   theme_bw() +
#   facet_grid(. ~ Group, scale = "free_x", space = "free_x") +
#   scale_y_log10(breaks = round(c(1/8, 1/4, 1/2, 1, 2, 4, 8), 2)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#   plot.title = element_text(hjust = 0.5)) +
#   labs(x = "IPSSM",
#     color = "Sub-group") +
#   ggtitle("AMLt (Reference: IPSSM Moderate-High)")

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
# hr_aml_group_plot
# dev.off()

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
  ggplot(aes(y = IPSSM, x = P*100, fill = sub_group, color = sub_group)) +
  geom_point(stat = "identity", position = "dodge") +
  geom_errorbar(aes(xmin = CILOW*100, xmax = CIHIGH*100)) +
  facet_grid2(sub_group ~ Time, scales = "free_y", 
  switch = "y", 
  space = "free_y",
    strip = strip_themed(
                        background_y = elem_list_rect(fill = colors),
                         text_y = element_text(angle = 0, color = "white", face = "bold")
                )) +
  scale_fill_manual(values = c(colors, "black")) +
  scale_color_manual(values =  c(colors, "black")) +
  scale_y_discrete(name = "", position = "right") +
  theme_bw() +
  theme( plot.title = element_text(hjust = 0.5), legend.position = "none",
  axis.x.text = element_text(angle = 90, hjust = 1)) +
  labs(x = "", x = "Probability of AMLt",
    fill = "", color = "") +
  ggtitle("AMLt probability")


png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_P_subgroups2.png", width = 2000, height = 2500, res = 300)
AMLt_P_plot
dev.off()


## Compute effect of AML
coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX,  data = IWS_full_amlt)

# getCoefsAML <- function(df, group, IPSSM_cats = levels(df$IPSSM)){
#   df_sub <- df %>%
#     filter(sub_group == group & IPSSM %in% IPSSM_cats) 
#   mod <- summary(coxph(Surv(AMLt_YEARS,PROG_STATE == "AMLt") ~ IPSSM_SCORE + AGE + SEX, data = df_sub))
#   mod$coefficients[1, c(1, 5)]
# }
# ipssm_amlt_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt)
# ipssm_amlt_lowRisk_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt, 
#   IPSSM = c("Very-Low", "Low", "Moderate-Low"))
# ipssm_amlt_HighRisk_effect_list <- lapply(groups, getCoefsAML, df = IWS_amlt, 
#   IPSSM = c("Moderate-High", "High", "Very-High"))

# ipssm_amlt_effect_list_joint <- c(ipssm_amlt_effect_list, ipssm_amlt_lowRisk_effect_list, ipssm_amlt_HighRisk_effect_list)

# ipssm_amlt_effect_tib <- tibble(group = names(ipssm_amlt_effect_list_joint), 
#   IPSSM_effect = sapply(ipssm_amlt_effect_list_joint, function(x) x[1]),
#   p_value = sapply(ipssm_amlt_effect_list_joint, function(x) x[2]),
#   Prog = "AMLt", Category = rep(c("All", "Low-risk", "High-risk"), each = length(groups)))

# ipssm_prog_joint <- bind_rows(ipssm_effect_tib, ipssm_amlt_effect_tib)
# write.table(ipssm_prog_joint, file = "results/GESMD_IWS_clustering/IPSSM_effects_OS_AMLt.txt", 
#   sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

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


# ## TET2-bi effect in low-risk IPSSM
# amlt_TET2bi_morph <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ TET2bi_group, 
#   IWS_amlt, subset = IPSSM %in% c("Very-Low", "Low", "Moderate-Low")) %>%
#     ggsurvplot(data = IWS_amlt,  palette = c("black", "#56B4E9"),
#      risk.table = TRUE, break.time.by = 2, fun = "event", xlim = c(0, 6),
#      ylim = c(0, 0.2),
#       legend.labs  = c("Other", "TET2-bi"))

# aml_tet2bi_morph <- plot_grid(amlt_TET2bi_morph$plot + 
#         ggtitle("IWS Morphologic Low Risk") +
#          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#          ylab("AMLt probability") +
#          xlab("Time (years)"), 
#         amlt_TET2bi_morph$table, ncol = 1)
# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_tet2bi_morph_lowrisk.png", width = 2000, height = 2000, res = 300)
# aml_tet2bi_morph
# dev.off()

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


# ## AMLt by IPSSM (only IWS)
# amlt_ipssm <- lapply(IPSSM_groups, function(cat){
#   df <- filter(IWS_mds_f, IPSSM == cat)
#   n_group <- df %>% 
#     group_by(sub_group) %>%
#     summarize(n = n())
#   sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
#   df <- filter(df, sub_group %in% sel_clusts)
#   p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group, df) %>%
#     ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
#                palette = colors[which(groups %in% sel_clusts)],
#                risk.table = TRUE, break.time.by = 2, xlim = c(0, 10), 
#                xlab = "Time (Years)")
#   p
# })
# amlt_ipssm_plots <- lapply(IPSSM_groups, function(ipssm){
#    aml_plot <-  plot_grid(
#         plot_grid(amlt_ipssm[[ipssm]]$plot + 
#             ggtitle(ipssm) +
#              theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#              ylab("AMLt probability") +
#              xlab("Time (years)") +
#              ylim(0, 1), 
#             amlt_ipssm[[ipssm]]$table, ncol = 1),
#         ncol = 1)
#     # ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_", ipssm, ".png"),
#     #     width = 2000, height = 2000, dpi = 300, units = "px")
#     aml_plot
# })

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/AMLt_subgroups_IPSSM_panel.png", width = 4000, height = 6000, res = 300)
# plot_grid(plotlist = amlt_ipssm_plots, ncol = 2, labels = "AUTO")
# dev.off()

# ## AMlt by subgroup (only IWS)
# amlt_subgroups <- lapply(groups, function(group){
#   df <- filter(IWS_mds, sub_group == group)
#   n_ipssm <- df %>% 
#     group_by(IPSSM) %>%
#     summarize(n = n())
#   sel_ipssm <- as.character(filter(n_ipssm, n >= 10)$IPSSM)
#   df <- filter(df, IPSSM %in% sel_ipssm)
#   p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ IPSSM, df) %>%
#     ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
#                palette = ipssm_cols[which(IPSSM_groups %in% sel_ipssm)],
#                risk.table = TRUE, break.time.by = 2) +
#     xlab("Time (Years)")
#   p
# })
# names(amlt_subgroups) <- groups
# amlt_subgroups_plots <- lapply(names(amlt_subgroups), function(group){
#     aml_plot <- plot_grid(
#         plot_grid(amlt_subgroups[[group]]$plot + 
#             ggtitle(group) +
#              theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#              ylab("AMLt probability") +
#              xlab("Time (years)") +
#              ylim(0, 1), 
#             amlt_subgroups[[group]]$table, ncol = 1),
#         ncol = 1)
#     png_filename <- paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_", group, ".png")
#     png_filename <- gsub("[^A-Za-z0-9/_-].", "_", png_filename, fixed = TRUE)
#     ggsave(png_filename, width = 2000, height = 2000, dpi = 300, units = "px")
#     aml_plot
# })
# png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_subgroups_panel.png", width = 4000, height = 6000, res = 300)
# plot_grid(plotlist = amlt_subgroups_plots, ncol = 2, labels = "AUTO")
# dev.off()


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
  Gene2 = c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk", "Cytogenetic Risk", "DNMT3A"), 
  x = c("EZH2", "STAG2", "EZH2", "STAG2", "EZH2", "STAG2", "TET2-bi"),           
  y = c(4.4),                      
  label = c("**", "*", "*", "**", "*", "**", "*")  
)



os_int_plot_main <- mut_hr_comb_joint %>%
  filter(!is.na(Gene) & Gene %in% c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk", "DNMT3A")) %>%
  mutate(Gene2 = fct_relevel(droplevels(Gene), c("Cytogenetic Risk", "+8", "U2AF1", "SETBP1", "TET2 mono-allelic", "DNMT3A"))) %>%
  ggplot(aes(x = sub_group, y = HR, color = sub_group, fill = sub_group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), color = "black", width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-group") +
  ggtitle("Overall Survival") +
  facet_wrap(~ Gene2, scales = "free_x") +
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
  filter(!is.na(Gene) & !Gene %in% c("TET2 mono-allelic", "U2AF1", "SETBP1", "+8", "Cytogenetic Risk", "BM_BLAST", "HB", "PLT", "DNMT3A", "IDH1", "TP53 mono-allelic", "del20q")) %>%
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
  mono7 = factor(mono7, levels = c("-7", "Other")),
  tet2bi_group = ifelse(ID %in% filter(joint_mds, sub_group == "TET2-bi")$ID, "TET2-bi", "Other"),
  tet2bi_group = factor(tet2bi_group, levels = c("TET2-bi", "Other"))

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
summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ tet2bi_group*DNMT3A + AGE + SEX + dataset, joint_full_subs))

summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ stag2_group*RUNX1 + AGE + SEX + dataset, joint_full_subs)) ## No signif


summary(coxph(Surv(OS_YEARS, OS_STATUS) ~ as.numeric(CYTO_IPSSR) + AGE + SEX + dataset, subset = sub_group == "STAG2", data = joint_mds))

stag2_cyto <- 
    survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ CYTO_IPSSR, joint_mds, 
    subset = sub_group == "STAG2" & CYTO_IPSSR %in% c("Good", "Int")) %>%
        ggsurvplot(data = joint_mds, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, 
                 legend = "none", title = "Cytogenetic Risk in STAG2 subgroup",
                  pval = sprintf("HR = %.2f\nP = %.3f", 0.994970, 0.98201),
                 pval.coord = c(7, 0.9), legend.labs = c("Good", "Int"),
                 xlim = c(0, 10), palette = c("#66bd63", "#fee08b" ), 
        xlab = "Time (Years)") 

stag2_cyto$plot <- stag2_cyto$plot +  theme(plot.title = element_text(hjust = 0.5))

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_cyto.png", width = 2000, height = 2000, res = 300)
# stag2_cyto
# dev.off()

full_cyto <- 
    survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ CYTO_IPSSR, joint_full) %>%
        ggsurvplot(data = joint_full, surv.median.line = "hv",
                 risk.table = TRUE, break.time.by = 2, legend = "none",
                 title = "Cytogenetic Risk in full cohorts",
                 legend.labs = c("Very Good", "Good", "Int", "Poor", "Very Poor"),
                 xlim = c(0, 10), palette = c("#2ca25f", "#66bd63", "#fee08b", "#f46d43", "#d73027"), 
                 xlab = "Time (Years)") 
full_cyto$plot <- full_cyto$plot +  theme(plot.title = element_text(hjust = 0.5))

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/full_cyto.png", width = 2000, height = 2000, res = 300)
# full_cyto
# dev.off()



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

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/forest_cyto_stag2.png", width = 2000, height = 2000, res = 300)
# print(stag2_cyto_forest)
# dev.off()

full_cyto_tab <- make_forest_table(var_name = "as.numeric(CYTO_IPSSR)", data_list = data_list_full, 
  col1_fun = function(df) sum(!is.na(df$CYTO_IPSSR)), 
  col2_fun = function(df) return(NA),
  headers = c("N CYTO", "A"))

full_cyto_forest <- plot_forest(full_cyto_tab[, -3], plot_title = "Cytogenetic Risk in full dataset",
   summary_rows = summary_rows, ci_column = 3)

# png("figures/GESMD_IWS_clustering/subgroup_prog_inter/forest_cyto_full.png", width = 2000, height = 2000, res = 300)
# print(full_cyto_forest)
# dev.off()



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

## TET2-bi vs DNMT3A
data_list_tet2bi <- list(
    IWS = filter(joint_mds, dataset == "IWS" & sub_group == "TET2-bi"),
    GESMD = filter(joint_mds, dataset == "GESMD" & sub_group == "TET2-bi"),
    Joint = filter(joint_mds, sub_group == "TET2-bi"),
    Hershberger = filter(hersh_mds_cyto, sub_group == "TET2-bi") %>% mutate(TET2other = TET2mono)
)
tet2_dnmt3a_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("TET2-bi")), variable = "DNMT3A", title = "DNMT3A in TET2-bi subgroup")
tet2_dnmt3a_tab <- make_forest_table(var_name = "DNMT3A", data_list = data_list_tet2bi, 
  col1_fun = function(df) sum(df$DNMT3A == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$DNMT3A == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

tet2_dnmt3a_forest <- plot_forest(tet2_dnmt3a_tab[1:4, ], plot_title = "DNMT3A in TET2-bi sub-group",
   summary_rows = summary_rows[1:4])

full_dnmt3a_surv <- makeSurvPlot(joint_full, variable = "DNMT3A", title = "DNMT3A in full cohort")
full_dnmt3a_tab <- make_forest_table(var_name = "DNMT3A", data_list = data_list_full, 
  col1_fun = function(df) sum(df$DNMT3A == 0, na.rm = TRUE), 
  col2_fun = function(df) sum(df$DNMT3A == 1, na.rm = TRUE), 
  headers = c("WT", "MUT"))

full_dnmt3a_forest <- plot_forest(full_dnmt3a_tab, plot_title = "DNMT3A in full cohort",
   summary_rows = summary_rows)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/TET2bi_DNMT3A_panel.png", width = 3500, height = 2500, res = 300)
plot_grid(tet2_dnmt3a_surv, full_dnmt3a_surv, tet2_dnmt3a_forest, full_dnmt3a_forest, 
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

## STAG2 vs BM_BLAST
colores_blasts <- c("#c9bb92", "#FFB300", "#D32F2F")
# data_list_mono7 <- list(
#     IWS = filter(joint_mds, dataset == "IWS" & sub_group == "-7"),
#     GESMD = filter(joint_mds, dataset == "GESMD" & sub_group == "-7"),
#     Joint = filter(joint_mds, sub_group == "-7"),
#     Hershberger = filter(hersh_mds, sub_group == "-7") 
# )

# mono7_blast_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("-7")) %>% mutate(BLAST = cut(BM_BLAST, breaks = c(0, 5, 10, 40))), 
#   variable = "BLAST", title = "BM BLAST in -7 subgroup", legend.labs = c("<5%", "5-10%", "\\>10%"), palette = colores_blasts)
# mono7_blast_tab <- make_forest_table(var_name = "BM_BLAST", data_list = data_list_mono7,
#   col1_fun = function(df) sum(!is.na(df$BM_BLAST)), 
#   col2_fun = function(df) NA, 
#   headers = c("N BM BLAST", "A"))
# mono7_blast_forest <- plot_forest(mono7_blast_tab[, -3], plot_title = "BM BLAST in -7 subgroup",
#    summary_rows = summary_rows, ci_column = 3)


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

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_blast_panel.png", width = 3500, height = 2500, res = 300)
# plot_grid(
#   plot_grid(stag2_blast_surv, mono7_blast_surv, stag2_blast_forest, mono7_blast_forest, 
#     ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D")),
#   plot_grid(full_blast_surv, full_blast_forest, ncol = 1, rel_heights = c(2, 1), labels = c("E", "F")),
#   ncol = 1, rel_heights = c(1, 1)
# )
plot_grid(stag2_blast_surv, full_blast_surv, stag2_blast_forest, full_blast_forest, 
  ncol = 2, rel_heights = c(1.5, 1), labels = c("A", "C", "B", "D"))
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
ezh2_cyto_surv <- makeSurvPlot(subset(joint_mds, sub_group %in% c("EZH2") & !CYTO_IPSSR %in% c("Very-Good", "Poor")), 
  variable = "CYTO_IPSSR", title = "Cytogenetic Risk in EZH2 subgroup", legend.labs = c("Good", "Int"), 
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
summary(coxph(Surv(AMLt_YEARS, PROG_STATE == "AMLt") ~ stag2_group*plus8 + AGE + SEX, IWS_aml_groups)) # No signif


amlt_int_sig <- data.frame(
  Gene = c("ASXL1", "TET2 mono-allelic", "RUNX1", "ASXL1", "RUNX1", "BCOR", "+8"), 
  x = c("EZH2", "EZH2", "EZH2", "STAG2", "STAG2", "STAG2", "STAG2"),           
  y = c(0.15),                      
  label = c("**", "**", "*", "*", "***", "**", "*")  
)

aml_int_plot_main <- mut_hr_aml_comb_joint %>%
  filter(!is.na(Gene) & Gene %in% c("TET2 mono-allelic", "ASXL1", "RUNX1", "BCOR", "+8")) %>%
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
  filter(!is.na(Gene) & !Gene %in% c("CBL", "IDH1", "TET2 mono-allelic", "ASXL1", "RUNX1", "BCOR", "-Y", "BM_BLAST", "HB", "PLT", "TP53 mono-allelic", "del20q", "+8")) %>%
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



png("figures/GESMD_IWS_clustering/subgroup_prog_inter/prognosis_panel.png", width = 3200, height = 5000, res = 300)
plot_grid(
  plot_grid(median_OS_plot, AMLt_P_plot, 
  ncol = 1, rel_heights = c(1, 1.2), labels = c("A", "E")),
  plot_grid(os_int_plot_main,
      plot_grid(stag2_cyto$plot, stag2_cyto$table, ncol = 1, rel_heights = c(1.6, 1)),
  stag2_cyto_forest, aml_tet2bi_full, aml_int_plot_main, ncol = 1,
  rel_heights = c(1.3, 1.5, 1, 1.5, 1.3), labels = c("B", "C","D", "F", "G")),
  ncol = 2, rel_widths = c(1, 1.3)
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

## STAG2 and +8
stag2_plus8_aml <- makeSurvPlotAML(subset(IWS_amlt, sub_group == "STAG2"), 
  variable = "plus8", title = "+8 in STAG2 sub-group", pval_height = 0.10)

full_plus8_aml <- makeSurvPlotAML(IWS_amlt, 
  variable = "plus8", title = "+8 in IWS cohort", pval_height = 0.10)

png("figures/GESMD_IWS_clustering/subgroup_prog_inter/STAG2_plus8_aml_panel.png", res = 300, height = 1300, width = 2600)
plot_grid(stag2_plus8_aml, full_plus8_aml, ncol = 2, labels = "AUTO")
dev.off()

