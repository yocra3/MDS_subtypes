#' ---------------------------
#'
#' Purpose of script:
#'
#'  Make plots for EHA abstract
#' 
#' ---------------------------
#'
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(MASS)
library(tidyverse)
library(cowplot)
library(survminer)
library(survival)

load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

## Clinical summary plot
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_joint <- bind_rows(
    IWS_mds %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "GESMD")
) %>%
    mutate(sub_group = fct_recode(sub_group, 
    "del5q-IB" = "del5q",
    "SF3B1-IB" = "SF3B1",
    "TET2-bi" = "TET2 bi-allelic",
    "MDS-LB" = "Low blasts",
    "-7" = "7-"))

clin_sum_plot <- clin_joint %>% 
    select(all_of(clin_vars), sub_group, dataset) %>%
    pivot_longer(cols = all_of(clin_vars), names_to = "Clin_var", values_to = "Value") %>%
    filter(Clin_var != "PB_BLAST") %>%
    group_by(Clin_var) %>%
    mutate(sd_total = sd(Value, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(Clin_var, sub_group) %>%
    summarize(m = mean(Value, na.rm = TRUE),
        sd_total = first(sd_total)) %>% 
    ungroup() %>%
    group_by(Clin_var) %>%
    mutate(d = (m - m[ sub_group == "MDS-LB"])/sd_total,
    Clin_var = factor(Clin_var, levels = c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT"))) %>%
    filter(sub_group != "MDS-LB") %>%
    ggplot(aes(x = sub_group, y = Clin_var, fill = d )) +
    geom_tile() +
    scale_fill_gradient2(
    low = "#0033ff",  
    high = "#ff0000",
    midpoint = 0
    ) +
  ylab("Clinical") +
  xlab("Sub-groups") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1))

png("figures/abstract_EHA_2026/clinical_summary.png", width = 1500, height = 800, res = 300)
clin_sum_plot
dev.off()

## Mutations
mutations <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "CBL", "IDH1")

sub_groups <- c("EZH2", "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex", "MDS-LB", "MDS-IB1", "MDS-IB2")


joint_mut <- rbind(
    IWS_mds %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")))  %>%
    mutate(sub_group = fct_recode(sub_group, 
    "del5q-IB" = "del5q",
    "SF3B1-IB" = "SF3B1",
    "TET2-bi" = "TET2 bi-allelic",
    "MDS-LB" = "Low blasts",
    "-7" = "7-"))

mut_ORs <- lapply(mutations, function(var){
    lapply(sub_groups[sub_groups != "MDS-LB"], function(group){
        sel_df <- filter(joint_mut, sub_group %in% c(group, "MDS-LB")) %>%
            mutate(sub_group = factor(sub_group, levels = c("MDS-LB", group)))
        tab <- table(sel_df[[var]], sel_df$sub_group)
        fisher_res <- fisher.test(tab)
        tibble(Gene = var,
               Comp_clust = group,
               OR = fisher_res$estimate,
               P_value = fisher_res$p.value)
    }) %>%    Reduce(., f = rbind)
}) %>% Reduce(., f = rbind)



mut_sum_plot <- mut_ORs %>%
    filter(Comp_clust != "MDS-LB") %>%
    mutate(Gene = ifelse(Gene == "TET2other", "TET2-mono", Gene),
        Mut = factor(Gene, levels = c("RUNX1", "ASXL1", "ZRSR2", "SETBP1", "CBL", 
    "IDH2", "BCOR", "TET2-mono", "U2AF1", "IDH1", "DNMT3A", "SRSF2", "DDX41")), 
    sub_group = factor(Comp_clust, levels = sub_groups)) %>%
    mutate(OR = ifelse(OR == 0, 0.09, OR)) %>%
    ggplot(aes(x = sub_group, y = Mut, fill = log(OR) )) +
    geom_tile() +
    scale_fill_gradient2(
    low = "#0033ff",  
    mid = "white",    
    high = "#ff0000",  
    midpoint = 0,
    name = "",
    breaks = c(-2, 0, 2),
    labels = c("Group<LB", "Group=LB", "Group>LB")
    ) +
  ylab("Mutations") +
  xlab("Sub-groups") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1))

png("figures/abstract_EHA_2026/mutation_subgroups_summary.png", width = 1600, height = 1200, res = 300)
mut_sum_plot
dev.off()


## Karyotypes
kar_events <- list("plus8" = "+8", "delY" = "-Y", "del20q" = "del20q")

joint_kar <- rbind(
    IWS_mds %>% select(all_of(names(kar_events)), sub_group, CYTO_IPSSR) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(names(kar_events)), sub_group, CYTO_IPSSR) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
    CYTO_IPSSR = fct_collapse(CYTO_IPSSR,
    "Very-Good" = c("Very-Good", "Very Good"),
    "Intermediate" = c("Intermediate", "Int"),
    "Very-Poor" = c("Very-Poor", "Very Poor"))) %>%
    mutate(sub_group = fct_recode(sub_group, 
    "del5q-IB" = "del5q",
    "SF3B1-IB" = "SF3B1",
    "TET2-bi" = "TET2 bi-allelic",
    "MDS-LB" = "Low blasts",
    "-7" = "7-"))

kar_ORs <- lapply(names(kar_events), function(var){
    lapply(sub_groups[sub_groups != "MDS-LB"], function(group){
        sel_df <- filter(joint_kar, sub_group %in% c(group, "MDS-LB")) %>%
            mutate(sub_group = factor(sub_group, levels = c("MDS-LB", group)))
        tab <- table(sel_df[[var]], sel_df$sub_group)
        fisher_res <- fisher.test(tab)
        tibble(Event = kar_events[[var]],
               Sub_group = group,
               OR = fisher_res$estimate,
               P_value = fisher_res$p.value)
    }) %>%    Reduce(., f = rbind)
}) %>% Reduce(., f = rbind)


kar_sum_plot <- kar_ORs %>%
    filter(Sub_group != "MDS-LB") %>%
    mutate(Event = factor(Event, levels = c("+8", "-Y", "del20q")), 
    sub_group = factor(Sub_group, levels = sub_groups),
    OR = ifelse(OR == 0, 0.15, OR)) %>%
    ggplot(aes(x = sub_group, y = Event, fill = log(OR) )) +
    geom_tile() +
    scale_fill_gradient2(
    low = "#0033ff",   
    high = "#ff0000", 
    midpoint = 0
    ) +
  ylab("Karyotype") +
  xlab("Sub-groups") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1))

png("figures/abstract_EHA_2026/karyotype_subgroups_summary.png", width = 1400, height = 600, res = 300)
kar_sum_plot
dev.off()


## IPSSM

ipssm_summary <- rbind(
    IWS_mds %>% select(IPSSM, sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(IPSSM, sub_group) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
    mutate(sub_group = fct_recode(sub_group, 
    "del5q-IB" = "del5q",
    "SF3B1-IB" = "SF3B1",
    "TET2-bi" = "TET2 bi-allelic",
    "MDS-LB" = "Low blasts",
    "-7" = "7-")) %>%
  filter(!is.na(IPSSM)) %>%
  filter(sub_group %in% c("EZH2", "TET2-bi", "-7", "STAG2")) %>%
  group_by(sub_group, IPSSM, dataset) %>%
  summarize(N = n()) %>%
  group_by(sub_group, dataset) %>%
  mutate(Freq = N/sum(N),
         IPSSM = fct_recode(IPSSM, 
         "VL" = "Very-Low",
         "L" = "Low",
         "ML" = "Moderate-Low", 
         "MH" = "Moderate-High",
         "H" = "High",
         "VH" = "Very-High")) %>%
  ggplot(aes(x = dataset, y = Freq*100, fill = fct_rev(IPSSM))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "IPSSM classification",
    y = "Individuals (%)",
    x = "Sub-group",
    fill = "IPSSM") +
    facet_grid(~ sub_group) +
  theme(plot.title = element_text(hjust = 0.5), 
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) 

png("figures/abstract_EHA_2026/IPSSM_summary.png", res = 300, height = 700, width = 1400)
ipssm_summary
dev.off()



## OS STAG2

joint_mds <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS"),
    gesmd_dataset %>% mutate(dataset = "GESMD") 
)

joint_prognosis_plot <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS") ,
    gesmd_dataset %>% mutate(dataset = "GESMD")) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
     IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 



ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")

prognosis_stag2 <- joint_prognosis_plot %>% 
    filter(sub_group == "STAG2") %>%
    mutate(IPSSM_short = fct_recode(IPSSM, 
        "VL" = "Very-Low", "L" = "Low", "ML" = "Moderate-Low", 
        "MH" = "Moderate-High", "H" = "High", "VH" = "Very-High"))

stag2_plot <- ggsurvplot(survfit(Surv(OS_YEARS, OS_STATUS) ~ IPSSM_short, data = prognosis_stag2), 
    surv.median.line = "hv", palette = ipssm_cols,
        risk.table = TRUE, break.time.by = 2, conf.int = FALSE,
        xlim = c(0, 6), ylim = c(0, 1),
        legend.labs  = levels(prognosis_stag2$IPSSM_short),
        legend = "right", 
        legend.title = "IPSS-M")

png("figures/abstract_EHA_2026/STAG2_IPSSM_OS.png", width = 1200, height = 700, res = 300)
stag2_plot$plot + 
    ggtitle("IPSSM in STAG2") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("OS prob.") +
    xlab("Time (years)") 
dev.off()


plus8_plot <- ggsurvplot(survfit(Surv(OS_YEARS, OS_STATUS) ~ plus8, data = prognosis_stag2), 
    surv.median.line = "hv", palette = c("Blue", "Red"),
        risk.table = TRUE, break.time.by = 2, conf.int = FALSE,
        xlim = c(0, 6), ylim = c(0, 1),
        legend.labs  = c("WT", "+8"),
        legend = "right", 
        legend.title = "+8")

png("figures/abstract_EHA_2026/STAG2_plus8_OS.png", width = 1200, height = 700, res = 300)
plus8_plot$plot + 
    ggtitle("+8 in STAG2") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("OS prob.") +
    xlab("Time (years)")
dev.off()


prognosis_bitet2 <- IWS_mds %>% 
    filter(sub_group == "TET2 bi-allelic") %>%
    mutate(IPSSM_short = fct_recode(IPSSM, 
        "VL" = "Very-Low", "L" = "Low", "ML" = "Moderate-Low", 
        "MH" = "Moderate-High", "H" = "High", "VH" = "Very-High"))


tet2bi <- ggsurvplot(survfit(Surv(AMLt_YEARS, AMLt_STATUS) ~ IPSSM_short, data = prognosis_bitet2), 
        palette = ipssm_cols,
        fun = "event", break.time.by = 2, conf.int = FALSE,
        xlim = c(0, 6), ylim = c(0, 1),
        legend.labs  = levels(prognosis_bitet2$IPSSM_short),
        legend = "right", 
        legend.title = "IPSS-M")

png("figures/abstract_EHA_2026/TET2bi_IPSSM_OS.png", width = 1200, height = 700, res = 300)
tet2bi$plot + 
    ggtitle("IPSSM in TET2-bi") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("AMLt cum. inc.") +
    xlab("Time (years)")
dev.off()

all_IWS <- ggsurvplot(survfit(Surv(AMLt_YEARS, AMLt_STATUS) ~ IPSSM, data = IWS_mds), 
        palette = ipssm_cols,
        fun = "event", break.time.by = 2, conf.int = FALSE,
        xlim = c(0, 6), ylim = c(0, 1),
        legend.labs  = levels(IWS_mds$IPSSM),
        legend = "right", 
        legend.title = "IPSS-M")

png("figures/abstract_EHA_2026/all_IWS_IPSSM_OS.png", width = 1200, height = 700, res = 300)
all_IWS$plot + 
    ggtitle("IPSSM in all IWS") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("LFS prob.") +
    xlab("Time (years)")
dev.off()


## Clinical variables

## Clinical variables
clusters <- levels(IWS_mds$sub_group)


colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
    "#D55E00", "#999999", "grey40",  "black")
scale_col <- scale_color_manual(values = colors, name = "Sub-group")
scale_fill <- scale_fill_manual(values = colors, name = "Sub-group")

makeBoxPlots <- function(variable, label){
  rbind(gesmd_dataset %>% select(sub_group, !!sym(variable)) %>% mutate(dataset = "GESMD"),
      IWS_mds %>% select(sub_group, !!sym(variable)) %>% mutate(dataset = "IWS")) %>%
      mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
      ggplot(aes(x = dataset, y = !!sym(variable), fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab(label) +
      xlab("") +
      ggtitle(paste(label)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
}

png("figures/abstract_EHA_2026/clin_vars_boxplot.png", width = 1500, height = 800, res = 300)
clin_joint  %>%
    select(WBC, HB, PLT, sub_group) %>%
    mutate(PLT = ifelse(PLT > 423, 423, PLT),
           WBC = ifelse(WBC > 12, 12, WBC)) %>%
    pivot_longer(cols = c(WBC, HB, PLT), names_to = "Clin_var", values_to = "Value") %>%
    mutate(Clin_var = factor(Clin_var, levels = c("WBC", "HB", "PLT"))) %>%
    ggplot(aes(x = sub_group, y = Value, fill = sub_group)) +
    geom_boxplot() +
    facet_wrap(~ Clin_var, scales = "free_y") +
    scale_fill +
    ylab("Value") +
    xlab("") +
    ggtitle("Clinical variables") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
dev.off()


### Panel of clinical variables
clin_names <- list(WBC = "WBC count",
                   PLT = "Platelets",
                   HB = "Hemoglobin")
clin_plots <- lapply(names(clin_names), function(var){
    makeBoxPlots(var, clin_names[[var]]) +
    theme(legend.position = "none")
})

legend <- get_plot_component(makeBoxPlots("ANC", "ANC") + 
theme(legend.position = "bottom"),
 "guide-box", return_all = TRUE)[[3]]




png("figures/abstract_EHA_2026/clin_vars_panel.png", width = 1500, height = 800, res = 300)
plot_grid(
    plot_grid(plotlist = clin_plots, ncol = 3,
        labels = "AUTO"),
    legend,
    ncol = 1, rel_heights = c(1, 0.1))
dev.off()
