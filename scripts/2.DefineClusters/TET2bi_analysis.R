#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define prognosis for TET2 bi-allelic mutated patients
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
library(forestploter)

## Load data
load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/preprocess/clinical_preproc.Rdata")


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


plotFitTable <- function(df, mod_var, title, palette, labels,ylab, xlim ){

    mod <- formula(paste("Surv(AMLt_YEARS, AMLt_STATUS) ~", mod_var))
    my_fit <- survfit(formula = mod, df) 
    my_fit$call$formula <- mod
    plot <- ggsurvplot(my_fit, data = df,  palette = palette,
        risk.table = TRUE, break.time.by = 2, conf.int = FALSE,
        fun = "event", xlim = xlim, ylim = c(0, 1),
        legend.labs  = labels)


    plot_grid(plot$plot + 
        ggtitle(title) +
         theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         ylab(ylab) +
         xlab("Time (years)"), 
        plot$table, ncol = 1)
}
createTitle <- function(title)  { 
    ggplot() + 
  labs(title = title) +
  theme_void() + 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
}
ipssm_cols <-  c("#2ca25f",  "#66bd63", "#fee08b", "#fdae61",  "#f46d43", "#d73027")
tet2bi_col <- c("#D1EAF8", "#A3D5F1", "#56B4E9", "#4187AF", "#2C5B76", "#18303E")


iws_tet2 <- IWS_mds %>%
    filter(sub_group %in% c("TET2 bi-allelic", "Low blasts")) %>% 
    mutate(sub_group = droplevels(sub_group),
        sub_group_short = fct_recode(sub_group,   "LB" = "Low blasts", "TET2 bi" = "TET2 bi-allelic"),
        IPSSM_short = fct_recode(IPSSM, 
        "VL" = "Very-Low", "L" = "Low", "ML" = "Moderate-Low", 
        "MH" = "Moderate-High", "H" = "High", "VH" = "Very-High"),
    interaction = paste(sub_group_short, IPSSM_short, sep = " - "),
    interaction = factor(interaction, levels = c(paste(rep(c("LB", "TET2 bi"), each = 6), 
        c("VL", "L", "ML", "MH", "H", "VH"), sep = " - ")))) 

tet2_ipssm <- plotFitTable(iws_tet2, "interaction", "", c(ipssm_cols, tet2bi_col), 
    xlim = c(0, 15), labels = levels(iws_tet2$interaction), ylab = "AMLt free survival")
png("figures/GESMD_IWS_clustering/tet2bi/OS_interaction.png", width = 4500, height = 2000, res = 300)
plot_grid(createTitle("TET2 bi-allelic - IPSSM interaction"), tet2_ipssm, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

