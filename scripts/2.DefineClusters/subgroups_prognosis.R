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

load( "results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

## General survival
colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")

gesmd_cols <- colors[c(1:4, 11)]

clinical_blasts <- clinical_blasts %>%
    mutate(sub_group = factor(sub_group,
        levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", 
            "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2")))
gesmd_filt <- gesmd %>%
    filter(sub_group %in% c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "STAG2")) %>%
    mutate(sub_group = factor(sub_group,
        levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic", "STAG2")))

surv_IWS <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, clinical_blasts) %>%
    ggsurvplot(data = clinical_blasts, surv.median.line = "hv", palette = colors,
     risk.table = TRUE, break.time.by = 2, 
      legend.labs  = levels(clinical_blasts$sub_group))

surv_gesmd <- survfit(formula = Surv(OS_YEARS, OS_STATUS) ~ sub_group, gesmd_filt) %>%
    ggsurvplot(data = gesmd_filt, surv.median.line = "hv", 
               palette = gesmd_cols, risk.table = TRUE, break.time.by = 2, 
      legend.labs  = levels(gesmd_filt$sub_group)) 

png("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_all_subgroups.png", width = 4000, height = 2000, res = 300)
plot_grid(
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
    ncol = 2, rel_widths = c(1, 0.6))
dev.off()

## Survival by IPSSM (only IWS)
IPSSM_groups <- levels(clinical_blasts$IPSSM)
names(IPSSM_groups) <- IPSSM_groups
groups <- levels(clinical_blasts$sub_group)
names(groups) <- groups
survs_ipssm <- lapply(IPSSM_groups, function(cat){
  df <- filter(clinical_blasts, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = colors[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})
lapply(IPSSM_groups, function(ipssm){
    plot_grid(
        plot_grid(survs_ipssm[[ipssm]]$plot + 
            ggtitle(ipssm) +
             theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
             ylab("OS probability") +
             xlab("Time (years)"), 
            survs_ipssm[[ipssm]]$table, ncol = 1),
        ncol = 1)
    ggsave(paste0("figures/GESMD_IWS_clustering/subgroup_prognosis/OS_subgroups_", ipssm, ".png"),
        width = 2000, height = 2000, dpi = 300, units = "px")
})

### Moderate-High differences
iws_modhigh_os <- clinical_blasts %>%
    filter(IPSSM == "Moderate-High") %>%
    filter(sub_group %in% c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", "TET2 bi-allelic")) %>%
    mutate(TET2 = ifelse(sub_group %in% c("TET2 monoallelic", "TET2 bi-allelic"), "TET2", "Leukopenic")) 
summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ TET2 + AGE, data = iws_modhigh_os))

## Survival by subgroup (only IWS)
ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")

survs_subgroups <- lapply(groups[!groups %in% c("del20q", "del7q",	"complex")], function(group){
  df <- filter(clinical_blasts, sub_group == group)
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

## AML transformation

amlt_IWS <- survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ sub_group, clinical_blasts) %>%
    ggsurvplot(data = clinical_blasts, surv.median.line = "hv", palette = colors,
     risk.table = TRUE, break.time.by = 2, fun = "event",
      legend.labs  = levels(clinical_blasts$sub_group))

png("figures/GESMD_IWS_clustering/subgroup_prognosis/AMLt_all_subgroups.png", width = 2000, height = 2000, res = 300)
plot_grid(amlt_IWS$plot + 
        ggtitle("IWS") +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab("AMLt probability") +
         xlab("Time (years)"), 
        amlt_IWS$table, ncol = 1)
dev.off()

## AMLt by IPSSM (only IWS)
amlt_ipssm <- lapply(IPSSM_groups, function(cat){
  df <- filter(clinical_blasts, IPSSM == cat)
  n_group <- df %>% 
    group_by(sub_group) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$sub_group)
  df <- filter(df, sub_group %in% sel_clusts)
  p <- survfit(formula = Surv(AMLt_YEARS, AMLt_STATUS) ~ sub_group, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv", fun = "event",
               palette = colors[which(groups %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    xlab("Time (Years)")
  p
})
lapply(IPSSM_groups, function(ipssm){
    plot_grid(
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
})

## AMlt by subgroup (only IWS)
amlt_subgroups <- lapply(groups[!groups %in% c("del20q", "del7q",	"complex")], function(group){
  df <- filter(clinical_blasts, sub_group == group)
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
lapply(names(amlt_subgroups), function(group){
    plot_grid(
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
})