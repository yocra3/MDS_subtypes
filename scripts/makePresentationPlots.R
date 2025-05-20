library(MASS)
library(cowplot)
library(pheatmap)
library(survminer)
library(survival)
library(tidyverse)

load("./results/clustering/lowblast_filter_pc_clustering.Rdata")
load("./results/preprocess/clinical_preproc.Rdata")
load("./results/clustering/lowblast_filter_pc_clustering_raw.Rdata")

kar_events <- c("delY", "del11q", "del20q", "del7q",  "plus8",  "del7")
sel_muts <- c("TET2", "ASXL1", "SRSF2", "DNMT3A", "RUNX1", 
              "STAG2", "U2AF1", "EZH2", "ZRSR2", "TET2bi", "TET2other")
clinical_low$clust <- fct_recode(clinical_low$clust, `7-` = "del7", 
                                 `Highly Leukopenic` = "High Leukopenic", `Mildly Leukopenic` = "Mild Leukopenic")
clinical_low$logPLT <- log10(clinical_low$PLT)

### General PC ####
group_cols <- c("#000000", "#999999", "#E69F00",  
                "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"   
)
df_comb <- data.frame(comb_pc$x[, 1:4])
colnames(df_comb) <- paste0("PC", 1:4)
df_comb$cluster <- clinical_low$clust
p1 <- ggplot(df_comb, aes(x = PC1, y = PC2, col = cluster)) +
  geom_point() +
  xlab("PC1 (16.2%)") +
  ylab("PC2 (10.0%)") +
  scale_color_manual(values = group_cols) + 
  theme_bw() +
  theme(legend.position = "None")
p2 <- ggplot(df_comb, aes(x = PC3, y = PC4, col = cluster)) +
  geom_point() +
  xlab("PC3 (10.0%)") +
  scale_color_manual(values = group_cols) + 
  ylab("PC4 (9.0%)") +
  theme_bw()  +
  theme(legend.position = "None")

png("figures/MDS_2025/PCA.png", res = 300, heigh = 1200, width = 2000)
plot_grid(p1, p2, ncol = 2)
dev.off()

## Bone Marrow Blasts ####
png("figures/MDS_2025/BM_Blasts.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  mutate(category = cut(BM_BLAST, breaks = c(-0.1, 1, 2, 3, 4, 5, 10, 20), 
                        labels = c("0", "1", "2", "3", "4", "5-10", "10+"))) %>%
  group_by(category, clust) %>%
  summarize(n_cat = n()) %>%
  group_by(clust) %>%
  mutate(prop = n_cat/sum(n_cat),
         category = factor(category, levels = as.character(4:0))) %>%
  ggplot(aes(x = clust, y = prop*100, fill = category)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "Bone Marrow Blasts",
    y = "Individuals (%)",
    x = "Prognostic Group",
    fill = "BM (%)") +
  scale_fill_manual(values = c("#7a0177", "#c51b8a", "#f768a1", "#fa9fb5", "#fde0dd")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

## White Blood Cell counts ####
png("figures/MDS_2025/WBC.png", res = 300, heigh = 1200, width = 2000)
ggplot(clinical_low, aes(x = clust, y = WBC)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Whole Blood Cells ",
    y = expression(paste("WBC (x", 10^3, "/μL)")),
    x = "Prognostic Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

## Hemoglobin ####
png("figures/MDS_2025/HB.png", res = 300, heigh = 1200, width = 2000)

ggplot(clinical_low, aes(x = clust, y = HB)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Hemoglobin",
    y = "HB (g/dL)",
    x = "Prognostic Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

## Karyotypes ####

clinical_low$complex <- factor(clinical_low$complex, levels = c("non-complex", "complex"))
kar_events2 <- c("del7", "plus8", "delY", "del20q", "complex")
names(kar_events2) <- kar_events2
kar_tabs <- lapply(kar_events2, function(x) table(clinical_low[[x]], clinical_low$clust))
kar_tabs_lowblasts <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100

png("figures/MDS_2025/Karyotypes.png", res = 300, heigh = 1200, width = 2000)
pheatmap(mat = kar_tabs_lowblasts, display_numbers = TRUE)
dev.off()

# TET2 ####
## TET2 bi allelic
png("figures/MDS_2025/TET2bi.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  group_by(clust) %>%
  summarize(freq = mean(TET2bi)) %>%
  ggplot(aes(x = clust, y = freq*100, color = clust, fill = clust)) +
  geom_bar(stat = "identity") +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_bw() +
  labs(
    title = "TET2 bi-allelic",
    y = "Individuals (%)",
    x = "Sub-Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 16),
        legend.position = "none"
  )
dev.off()


a <- mutate(clinical_low, cot = ifelse(TET2bi == 1, "TET2bi",
                                       ifelse(TET2other == 1 & plus8 == 0, "TET2 mono",
                                              ifelse(del20q == 1, "del20q",
                                                     ifelse(delY == 1, "Y-",
                                                            ifelse(WBC > 5.5 & ANC > 2.5, "Mild",
                                                                   ifelse(del7 == 1, "7-",
                                                                          ifelse(plus8 == 1, "8+", "HL"))))))))
## TET2 monoallelic
png("figures/MDS_2025/TET2mono.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  group_by(clust) %>%
  summarize(freq = mean(TET2other)) %>%
  ggplot(aes(x = clust, y = freq*100, color = clust, fill = clust)) +
  geom_bar(stat = "identity") +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_bw() +
  labs(
    title = "TET2 mono-allelic",
    y = "Individuals (%)",
    x = "Sub-Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 16),
        legend.position = "none"
  )
dev.off()

### ASXL1 ####
png("figures/MDS_2025/ASXL1.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  group_by(clust) %>%
  summarize(freq = mean(ASXL1)) %>%
  ggplot(aes(x = clust, y = freq*100)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "ASXL1 mutations",
    y = "Individuals with mutation (%)",
    x = "Prognostic Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylim(c(0, 50))
dev.off()

### SRSF2 ####
png("figures/MDS_2025/SRSF2.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  group_by(clust) %>%
  summarize(freq = mean(SRSF2)) %>%
  ggplot(aes(x = clust, y = freq*100)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "SRSF2 mutations",
    y = "Individuals with mutation (%)",
    x = "Prognostic Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylim(c(0, 50))
dev.off()


### DNMT3A ####
png("figures/MDS_2025/DNMT3A.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  group_by(clust) %>%
  summarize(freq = mean(DNMT3A)) %>%
  ggplot(aes(x = clust, y = freq*100)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "DNMT3A mutations",
    y = "Individuals with mutation (%)",
    x = "Prognostic Group") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylim(c(0, 50))
dev.off()


## IPSSM ####
clinical_IPSSM <- mutate(clinical_low, 
                         IPSSM_col = ifelse(IPSSM %in% c("High", "Very-High"), "High/Very-High", 
                                            ifelse(grepl("Moderate", IPSSM), "Moderate",
                                                   ifelse(IPSSM == "Low", "Low", "Very-Low"))),
                         IPSSM_col = factor(IPSSM_col, levels = c("Very-Low", "Low", "Moderate", "High/Very-High")))
IPSSM_cats <- levels(clinical_IPSSM$IPSSM_col)
names(IPSSM_cats) <- IPSSM_cats

png("figures/MDS_2025/IPSSM_6cat.png", res = 300, heigh = 1200, width = 2000)
clinical_low %>%
  filter(!is.na(IPSSM)) %>%
  group_by(clust, IPSSM) %>%
  summarize(N = n()) %>%
  group_by(clust) %>%
  mutate(Freq = N/sum(N),
         clust = factor(clust, levels = c("Y-", "Mildly Leukopenic", "Highly Leukopenic", "TET2 mono", "TET2 bi", "del20q", "8+", "7-" ))) %>%
  ggplot(aes(x = clust, y = Freq*100, fill = fct_rev(IPSSM))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "IPSSM classification",
    y = "Individuals (%)",
    x = "Prognostic group",
    fill = "IPSSM Category") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  dev.off()

png("figures/MDS_2025/IPSSM_4cat.png", res = 300, heigh = 1200, width = 2000)
clinical_IPSSM %>%
  filter(!is.na(IPSSM_col)) %>%
  group_by(clust, IPSSM_col) %>%
  summarize(N = n()) %>%
  group_by(clust) %>%
  mutate(Freq = N/sum(N)) %>%
  ggplot(aes(x = clust, y = Freq*100, fill = fct_rev(IPSSM_col))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "IPSSM collapsed classification",
    y = "Individuals (%)",
    x = "Prognostic group",
    fill = "IPSSM Category") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dev.off()


## Survival ####
clinical_low2 <- mutate(clinical_low, 
                        cl = factor(clust, labels = c("HL", "ML", "T2m", "T2b", "Y-", "8+", "7-", "del20q")))


surv_low <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ cl, clinical_low2) %>%
  ggsurvplot(data = clinical_low2, surv.median.line = "hv",
             palette = group_cols, risk.table = TRUE, break.time.by = 2) +
  xlab("Time (Years)") +
  ggtitle("Overall Survival") 

png("figures/MDS_2025/Survival.png", res = 300, heigh = 1200, width = 2000)
surv_low$plot +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("figures/MDS_2025/Survival_table.png", res = 300, heigh = 700, width = 2000)
surv_low$table +
  theme(text = element_text(size = 12))
dev.off()

os_Ch_TET2bHL_df <- filter(clinical_low, 
                           clust %in% c("8+", "7-", "del20q", "TET2 bi", "High Leukopenic")) %>%
  mutate(cl = ifelse(clust %in% c("8+", "7-"), "Chromosomal", "Other"))
os_Ch_TET2bHL <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ cl + AGE, os_Ch_TET2bHL_df))

os_TET2bHL_TET2mML_df <- filter(clinical_low, 
                                clust %in% c("TET2 mono", "Mild Leukopenic", "TET2 bi", "High Leukopenic", "Y-", "del20q")) %>%
  mutate(cl = ifelse(clust %in% c("TET2 mono", "Mild Leukopenic", "Y-"), "Low", "Other"))
os_TET2bHL_TET2mML <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ cl + AGE, os_TET2bHL_TET2mML_df))

## Survival + IPSSM ####
clusts <- levels(clinical_low$clust)

IPSSM_groups <- levels(clinical_low$IPSSM)
names(IPSSM_groups) <- IPSSM_groups
survs <- lapply(IPSSM_groups, function(cat){
  df <- filter(clinical_low2, IPSSM == cat)
  n_group <- df %>% 
    group_by(clust) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$clust)
  df <- filter(df, clust %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ cl, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = group_cols[which(clusts %in% sel_clusts)],
               risk.table = TRUE, break.time.by = 2) +
    ggtitle(paste("IPSSM:", cat))  +
    xlab("Time (Years)")
  p
})



png("figures/MDS_2025/Survival_IPSSMLow.png", res = 300, heigh = 1200, width = 2000)
survs[["Low"]]$plot
dev.off()

png("figures/MDS_2025/Survival_IPSSMLow_table.png", res = 300, heigh = 700, width = 2000)
survs[["Low"]]$table +
  theme(text = element_text(size = 12))
dev.off()



png("figures/MDS_2025/Survival_IPSSMModerate_low.png", res = 300, heigh = 1200, width = 2000)
survs[["Moderate-Low"]]$plot
dev.off()


png("figures/MDS_2025/Survival_IPSSMModerate_low_table.png", res = 300, heigh = 500, width = 2000)
survs[["Moderate-Low"]]$table +
  theme(text = element_text(size = 12))
dev.off()


os_TET2b_ML_df <- filter(clinical_low, 
                         clust %in% c("TET2 mono", "Mildly Leukopenic", "TET2 bi", "Highly Leukopenic") & IPSSM == "Moderate-Low") %>%
  mutate(cl = ifelse(clust %in% c("TET2 bi"), "TET2 bi", "Other"))
os_TET2b_ML <- summary(coxph(Surv(OS_YEARS,OS_STATUS) ~ cl + AGE, os_TET2b_ML_df))



main_clusts <- c("TET2 mono", "Mildly Leukopenic", "TET2 bi", "Highly Leukopenic")
names(main_clusts) <- main_clusts

df_hl <- filter(clinical_low2, cl == "HL") %>%
  mutate(IPSSM = factor(IPSSM, labels = c("VL", "L", "ML", "MH", "H", "VH")))
p_hl <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ IPSSM, df_hl) %>%
  ggsurvplot(data = df_hl, surv.median.line = "hv",
               palette = c("#2ca25f", "#66bd63", "#fee08b", "#fdae61", "#f46d43", "#d73027"),
               risk.table = TRUE, break.time.by = 2) +
  ggtitle("HL") +
  xlab("Time (Years)")

df_t2b <- filter(clinical_low2, cl == "T2b" & IPSSM != "Very-High") %>%
  mutate(IPSSM = factor(IPSSM, labels = c("VL", "L", "ML", "MH", "H")))
p_t2b <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ IPSSM, df_t2b) %>%
  ggsurvplot(data = df_t2b, surv.median.line = "hv",
             palette = c("#2ca25f", "#66bd63", "#fee08b", "#fdae61", "#f46d43", "#d73027"),
             risk.table = TRUE, break.time.by = 2) +
  ggtitle("t2b") +
  xlab("Time (Years)")

png("figures/MDS_2025/Survival_HL.png", res = 300, heigh = 1200, width = 2000)
p_hl$plot
dev.off()


png("figures/MDS_2025/Survival_HL_table.png", res = 300, heigh = 600, width = 2000)
p_hl$table +
  theme(text = element_text(size = 12))
dev.off()


png("figures/MDS_2025/Survival_T2b.png", res = 300, heigh = 1200, width = 2000)
p_t2b$plot
dev.off()


png("figures/MDS_2025/Survival_T2b_table.png", res = 300, heigh = 500, width = 2000)
p_t2b$table +
  theme(text = element_text(size = 12))
dev.off()


survs3 <- lapply(levels(clinical_low$IPSSM), function(cat){
  df <- filter(clinical_low, IPSSM == cat)
  n_group <- df %>% 
    group_by(clust) %>%
    summarize(n = n())
  sel_clusts <- as.character(filter(n_group, n >= 10)$clust)
  df <- filter(df, clust %in% sel_clusts)
  p <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, df) %>%
    ggsurvplot(data = df, surv.median.line = "hv",
               palette = group_cols[which(clusts %in% sel_clusts)]) +
    ggtitle(paste("IPSSM:", cat))  
  p$plot +
    theme(plot.title = element_text(hjust = 0.5))
})


df_TET2bML_RestMH <- filter(clinical_low2, 
                            (clust == "TET2 bi" & IPSSM == "Moderate-Low") |
                              (clust %in% c("TET2 mono", "Mildly Leukopenic", "Highly Leukopenic") & IPSSM == "Moderate-High")) %>%
  mutate(clu = ifelse(cl == "T2b", "T2b","Other"))

p_TET2bML_RestMH <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ cl, df_TET2bML_RestMH) %>%
  ggsurvplot(data = df_TET2bML_RestMH, surv.median.line = "hv",
             palette = group_cols[1:4],
             risk.table = TRUE, break.time.by = 2)  +
  xlab("Time (Years)")


coxph(Surv(OS_YEARS,OS_STATUS) ~ clu + AGE, df_TET2bML_RestMH)


df_TET2bML_RestH <- filter(clinical_low2, 
                            (clust == "TET2 bi" & IPSSM == "Moderate-Low") |
                              (clust %in% c("TET2 mono", "Mildly Leukopenic", "Highly Leukopenic") & IPSSM == "High"))  %>%
  mutate(clu = ifelse(cl == "T2b", "T2b","Other"))
coxph(Surv(OS_YEARS,OS_STATUS) ~ clu + AGE, df_TET2bML_RestH)


p_TET2bML_RestH <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ cl, df_TET2bML_RestH) %>%
  ggsurvplot(data = df_TET2bML_RestH, surv.median.line = "hv",
             palette = group_cols[1:4],
             risk.table = TRUE, break.time.by = 2)  +
  xlab("Time (Years)")


png("figures/MDS_2025/Survival_T2B_ML_Rest_MH.png", res = 300, heigh = 1200, width = 2000)
p_TET2bML_RestMH$plot
dev.off()


png("figures/MDS_2025/Survival_T2B_ML_Rest_MH_table.png", res = 300, heigh = 500, width = 2000)
p_TET2bML_RestMH$table +
  theme(text = element_text(size = 12))
dev.off()

png("figures/MDS_2025/Survival_T2B_ML_Rest_H.png", res = 300, heigh = 1200, width = 2000)
p_TET2bML_RestH$plot
dev.off()


png("figures/MDS_2025/Survival_T2B_ML_Rest_H_table.png", res = 300, heigh = 500, width = 2000)
p_TET2bML_RestH$table +
  theme(text = element_text(size = 12))
dev.off()



## Mutations interaction ####
test_muts <- c("ASXL1", "SRSF2")

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
mut_hr_clust <- lapply(test_muts, function(gene){
  lapply(clusts, function(cl){
    df <- subset(clinical_low, clust == cl)
    df_est <- getEstimates(gene, df) %>%
      mutate(Cluster = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_ipssm <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, clinical) %>%
      mutate(Cluster = "IPSSM cohort")
  }) %>% Reduce(f = rbind)

mut_hr_plot <- rbind(mut_hr_clust, mut_hr_ipssm) %>%
  mutate(Cluster = factor(Cluster, levels = c(clusts, "IPSSM cohort"))) %>%
  filter(N >= 10) 
  

png("figures/MDS_2025/HR_mutations.png", res = 300, heigh = 1200, width = 2500)
mut_hr_plot %>%
  ggplot(aes(x = Cluster, y = HR, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = HR_L, ymax = HR_H), width = 0.2) +
  theme_bw() +
  scale_y_continuous(transform = "log2") +
  xlab("Sub-Group") +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_manual(name = "", values = c(group_cols[c(1:4, 6:7)], "white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 16))
dev.off()

## SRSF2 interaction
tet2bi_df <- clinical_low %>%
  mutate(ASXL1 = ifelse(ASXL1 == 1, "mutated", "WT"),
         SRSF2 = ifelse(SRSF2 == 1, "mutated", "WT")) %>%
  filter(clust == "TET2 bi")


p_T2b_SRSF2 <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SRSF2, tet2bi_df) %>%
  ggsurvplot(data = tet2bi_df, surv.median.line = "hv",
             risk.table = TRUE, break.time.by = 2) +
  xlab("Time (Years)")
png("figures/MDS_2025/T2B_SRSF2.png", res = 300, heigh = 1200, width = 2000)
p_T2b_SRSF2$plot
dev.off()

png("figures/MDS_2025/T2B_SRSF2_table.png", res = 300, heigh = 400, width = 2000)
p_T2b_SRSF2$table
dev.off()


clinical_df <- clinical %>%
  mutate(ASXL1 = ifelse(ASXL1 == 1, "mutated", "WT"),
         SRSF2 = ifelse(SRSF2 == 1, "mutated", "WT"),
         OS_YEARS = ifelse(OS_YEARS > 10, 10, OS_YEARS)) 

p_all_SRSF2 <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SRSF2, clinical_df) %>%
  ggsurvplot(data = tet2bi_df, surv.median.line = "hv",
             risk.table = TRUE, break.time.by = 2) +
  xlab("Time (Years)")

png("figures/MDS_2025/all_SRSF2.png", res = 300, heigh = 1200, width = 2000)
p_all_SRSF2$plot
dev.off()

png("figures/MDS_2025/all_SRSF2_table.png", res = 300, heigh = 400, width = 2000)
p_all_SRSF2$table
dev.off()

