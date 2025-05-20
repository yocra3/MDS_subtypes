#' ---------------------------
#'
#' Purpose of script:
#'
#'  Explore individual C-index for IPSSM
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.5 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
library(cowplot)
library(survival)
library(pheatmap)
library(survminer)


load("results/preprocess/clinical_c_index.Rdata")
load("results/clustering/umap_mutations_clustering.Rdata")



## Compute C-index per subtype
clinical_clust$RISK_IPSSM <- -clinical_clust$IPSSM_SCORE
concordance(Surv(OS_YEARS, OS_STATUS) ~ RISK_IPSSM, data = clinical_clust)
# Call:
# concordance.formula(object = Surv(OS_YEARS, OS_STATUS) ~ RISK_IPSSM, 
#     data = clinical_clust)

# n=2596 (257 observations deleted due to missingness)
# Concordance= 0.7483 se= 0.007228
# concordant discordant     tied.x     tied.y    tied.xy 
#    1276943     428497       2941        478          0 

subtypes <- unique(clinical_clust$consensus)
names(subtypes) <- subtypes
lapply(subtypes, function(subtype) {
    concordance(Surv(OS_YEARS, OS_STATUS) ~ RISK_IPSSM, 
        data = clinical_clust, subset = consensus == subtype)

})


comb2 <- left_join(clinical_clust,
    select(clinical_c_index, ID, OV_C_INDEX, SUB_C_INDEX), by ="ID")

group_by(clinical_c_index, consensus) %>%
    summarise(mean_ov50 = mean(OV_C_INDEX < 0.5, na.rm = TRUE),
             sum_ov50 = sum(OV_C_INDEX < 0.5, na.rm = TRUE),
              mean_sub50 = mean(SUB_C_INDEX < 0.5, na.rm = TRUE),
              sum_sub50 = sum(SUB_C_INDEX < 0.5, na.rm = TRUE))

group_by(comb2, consensus, cluster) %>%
    summarise(mean_ov50 = mean(OV_C_INDEX < 0.5, na.rm = TRUE),
             sum_ov50 = sum(OV_C_INDEX < 0.5, na.rm = TRUE),
              mean_sub50 = mean(SUB_C_INDEX < 0.5, na.rm = TRUE),
              sum_sub50 = sum(SUB_C_INDEX < 0.5, na.rm = TRUE))


# Plot
ov_plot <- ggplot(comb2, aes(x = OV_C_INDEX)) +
    geom_histogram() +
    facet_grid(consensus ~.) +
    theme_bw()
sub_plot <- ggplot(comb2, aes(x = SUB_C_INDEX)) +
    geom_histogram() +
    facet_grid(consensus ~.) +
    theme_bw()
ov_score_plot <- ggplot(comb2, aes(x = OV_C_INDEX, y = IPSSM_SCORE)) +
    geom_point() +
    facet_grid(consensus ~.) +
    theme_bw()
SUB_score_plot <- ggplot(comb2, aes(x = SUB_C_INDEX, y = IPSSM_SCORE)) +
    geom_point() +
    facet_grid(consensus ~.) +
    theme_bw()  

png("figures/overall_c_index.png", width = 1000)
plot_grid(ov_plot, sub_plot, ov_score_plot, SUB_score_plot, ncol = 4,
rel_widths = c(1, 1, 1, 1))
dev.off()

ov_plot_clus <- ggplot(comb2, aes(y = OV_C_INDEX, x = cluster)) +
    geom_boxplot() +
    facet_grid(~ consensus, scales = "free_x") +
    theme_bw()
sub_plot_clus <- ggplot(comb2, aes(y = SUB_C_INDEX, x = cluster)) +
    geom_boxplot() +
    facet_grid(~ consensus, scales = "free_x") +
    theme_bw()
ov_score_plot_clus <- ggplot(comb2, aes(x = OV_C_INDEX, y = IPSSM_SCORE, color = cluster)) +
    geom_point() +
    facet_grid(~ consensus) +
    theme_bw()
SUB_score_plot_clus <- ggplot(comb2, aes(x = SUB_C_INDEX, y = IPSSM_SCORE, color = cluster)) +
    geom_point() +
    facet_grid(~consensus ) +
    theme_bw()

png("figures/overall_clust_c_index.png", width = 1000, height = 1000)
plot_grid(ov_plot_clus, sub_plot_clus, ov_score_plot_clus, SUB_score_plot_clus, nrow = 4,
rel_heights = c(1, 1, 1.5, 1.5))
dev.off()

## Divide low-blast samples in two groups based on C-Index individual
lowb_data <- subset(comb2, consensus == "Low blasts" & !is.na(OV_C_INDEX) & !is.na(SUB_C_INDEX))
lowb_cind_ov_clusters <- kmeans(lowb_data$OV_C_INDEX, centers = 2)$cluster
lowb_cind_sub_clusters <- kmeans(lowb_data$SUB_C_INDEX, centers = 2)$cluster

tapply(lowb_data$OV_C_INDEX, lowb_cind_ov_clusters,  summary)
tapply(lowb_data$SUB_C_INDEX, lowb_cind_sub_clusters,  summary)

lowb_data$OV_CIND_CLUSTER <- factor(lowb_cind_ov_clusters)
lowb_data$SUB_CIND_CLUSTER <- factor(lowb_cind_sub_clusters)
lowb_data$OV_SUB_CIND_prop <- lowb_data$OV_C_INDEX/lowb_data$SUB_C_INDEX 
lowb_data$RATIO_CIND_CLUSTER <- factor(as.numeric(lowb_data$OV_SUB_CIND_prop > 1) + 1, levels = c(1, 2))
lowb_data <- mutate(lowb_data, 
    RATIO_CIND_CLUSTER = factor(ifelse(SUB_C_INDEX < 1.15*OV_C_INDEX - 0.13, 1, 2), levels = c(1, 2)),
    CIND_clust = factor(ifelse(SUB_CIND_CLUSTER == 1 & RATIO_CIND_CLUSTER == 1, "LowSub-LowConc",
        ifelse(SUB_CIND_CLUSTER == 1 & RATIO_CIND_CLUSTER == 2, "LowSub-HighConc",
            ifelse(SUB_CIND_CLUSTER == 2 & RATIO_CIND_CLUSTER == 1, "HighSub-LowConc", "HighSub-HighConc")))))


lowblast_ov <- ggplot(lowb_data, aes(x = OV_C_INDEX, y = SUB_C_INDEX, color = OV_CIND_CLUSTER)) +
    geom_point() +
    theme_bw()

lowblast_sub <- ggplot(lowb_data, aes(x = OV_C_INDEX, y = SUB_C_INDEX, color = SUB_CIND_CLUSTER)) +
    geom_point() +
    theme_bw()

lowblast_prop1 <- ggplot(lowb_data, aes(x = OV_SUB_CIND_prop)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, 2) +
    theme_bw()

lowblast4 <- ggplot(lowb_data, aes(x = OV_C_INDEX, y = SUB_C_INDEX, color = CIND_clust)) +
    geom_point() +
    theme_bw()


png("figures/lowblast_c_index_comp.png", width = 1000)
plot_grid(lowblast_ov, lowblast_sub, lowblast_prop1,  lowblast4, ncol = 2)
dev.off()


table(lowb_data$cluster, lowb_data$CIND_clust)
chis <- chisq.test(table(lowb_data$cluster, lowb_data$CIND_clust))
chis
chis$expected

coxph(Surv(OS_YEARS,OS_STATUS) ~ cluster*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 


## Sex
#' ---------------------------
table(lowb_data$CIND_clust, lowb_data$SEX)
chisq.test(table(lowb_data$CIND_clust, lowb_data$SEX))

## Age
#' ---------------------------
tapply(lowb_data$AGE, lowb_data$CIND_clust, summary)

## IPSSM
#' ---------------------------
table(lowb_data$CIND_clust, lowb_data$IPSSM)

ipssm1 <- ggplot(lowb_data, aes(x = OV_C_INDEX, y = SUB_C_INDEX, color = IPSSM)) +
    geom_point() +
    theme_bw()
ipssm2 <- ggplot(lowb_data, aes(y = OV_C_INDEX, x = IPSSM)) +
    geom_boxplot() +
    theme_bw()
ipssm3 <- ggplot(lowb_data, aes(y = SUB_C_INDEX, x = IPSSM)) +
    geom_boxplot() +
    theme_bw()

png("figures/lowblast_c_index_ipsmm.png", width = 1500)
plot_grid(ipssm1, ipssm2, ipssm3, ncol = 3)
dev.off()

png("figures/lowblast_c_index_ipsmm_cor.png", width = 1000)
ggplot(lowb_data, aes(x = OV_C_INDEX, y = SUB_C_INDEX, color = IPSSM)) +
    geom_point() +
    facet_wrap(~IPSSM) +
    theme_bw()
dev.off()


png("figures/lowblast_c_index_ipsmm_distr.png", width = 1000)
ipssm_d1 <- ggplot(lowb_data, aes(x = OV_C_INDEX, fill = IPSSM)) +
    geom_histogram() +
    facet_wrap(~IPSSM) +
    theme_bw()
ipssm_d2 <- ggplot(lowb_data, aes(x = SUB_C_INDEX, fill = IPSSM)) +
    geom_histogram() +
    facet_wrap(~IPSSM) +
    theme_bw()
plot_grid(ipssm_d1, ipssm_d2, nrow = 2)

dev.off()


## Kariotipos
#' ---------------------------
kar_events <- c("delY", "del11q", "del5q", "del12p",
                "del20q", "del7q", "plus8", "plus19", "del7", "complex")
names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(lowb_data[[x]], lowb_data$CIND_clust))
kar_tabs$complex <- kar_tabs$complex[2:1, ]

png("figures/lowblast_c_index_kariotypes.png")
kar_tabs_prop2 <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = kar_tabs_prop2, display_numbers = TRUE)
dev.off()
coxph(Surv(OS_YEARS,OS_STATUS) ~ plus8*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ del7*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ delY*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ complex*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 


## Mutaciones
#' ---------------------------
sel_muts <- c("TET2", "ASXL1", "DNMT3A", "RUNX1", "SF3B1", "SRSF2",  "U2AF1", "ZRSR2", "EZH2", "SETBP1", "STAG2")
names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(lowb_data[[x]], lowb_data$CIND_clust))

png("figures/lowblast_c_index_mutations.png", width = 700, height = 500)
mut_tabs_prop2 <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = mut_tabs_prop2, display_numbers = TRUE)
dev.off()

coxph(Surv(OS_YEARS,OS_STATUS) ~ TET2*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ STAG2*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*CIND_clust + AGE + SEX + IPSSR_SCORE, data = lowb_data) 


## Treatment
table(lowb_data$hma, lowb_data$IPSSM, lowb_data$CIND_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*CIND_clust + SEX + AGE + IPSSM_SCORE, data = lowb_data) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*IPSSM_SCORE + CIND_clust + SEX + AGE, data = lowb_data) 


cclust_low <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ CIND_clust + hma, lowb_data) %>%
    ggsurvplot(data = lowb_data)
png("figures/lowblast_c_index_surv.png", width = 1000, height = 500)
cclust_low
dev.off()