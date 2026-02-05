#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compare clustering subgroups with taxonomy classification
#' 
#' ---------------------------
#'
#' Notes:
#' Compare the clustering subgroups defined in the previous scripts with the
#' taxonomy classification of MDS patients. 
#' 
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------



# Load libraries and data
library(cowplot)
library(survminer)
library(survival)
library(tidyverse)
library(ipssm)
library(pheatmap)

load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/clustering/MFA_results.Rdata")

## Overlap with taxonomy
tax_tab <- table(IWS_mds$MOLECULAR_GROUP, IWS_mds$sub_group)

p_iws_tab <- pheatmap(mat = tax_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0, main = "IWS")


png("figures/GESMD_IWS_clustering/taxonomy_comparison/IWS_Taxonomy_mapping.png", width = 1800, height = 1800, res = 300)
p_iws_tab
dev.off()

tax_tab_gesmd <- table(gesmd_dataset$mol_manual, gesmd_dataset$sub_group)

p_gesmd_tab <- pheatmap(mat = tax_tab_gesmd, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0, main = "GESMD")

png("figures/GESMD_IWS_clustering/taxonomy_comparison/GESMD_Taxonomy_mapping.png", width = 1800, height = 1800, res = 300)
p_gesmd_tab
dev.off()


png("figures/GESMD_IWS_clustering/taxonomy_comparison/joint_Taxonomy_mapping.png", width = 3500, height = 1800, res = 300)
plot_grid(p_iws_tab$gtable, p_gesmd_tab$gtable , ncol = 2)
dev.off()

## Match manual clusters with original clusters
IWS_mfa_data <- left_join(IWS_dataset, 
    select(IWS_mds, ID, MOLECULAR_GROUP), by = "ID") 

IWS_tab <- table(IWS_mfa_data$MOLECULAR_GROUP, clust_comb1)

p_iws_clust <- pheatmap(mat = IWS_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0, main = "IWS")

png("figures/GESMD_IWS_clustering/taxonomy_comparison/IWS_cluster_mapping.png", width = 1800, height = 1500, res = 300)
p_iws_clust
dev.off()

GESMD_mfa_data <- left_join(gesmd_dataset_filt, 
    select(gesmd_dataset, ID, mol_manual), by = "ID") 


GESMD_tab <- table(GESMD_mfa_data$mol_manual, clust_gesmd)

p_gesmd_clust <- pheatmap(mat = GESMD_tab, display_numbers = TRUE, cluster_rows = FALSE, 
    cluster_cols = FALSE, legend = FALSE, angle_col = 0, main = "GESMD")

png("figures/GESMD_IWS_clustering/taxonomy_comparison/GESMD_cluster_mapping.png", width = 1800, height = 1500, res = 300)
p_gesmd_clust
dev.off()

png("figures/GESMD_IWS_clustering/taxonomy_comparison/joint_cluster_mapping.png", width = 3500, height = 1800, res = 300)
plot_grid(p_iws_clust$gtable, p_gesmd_clust$gtable , ncol = 2)
dev.off()


### IPSSM summary
ipssm_summary <- rbind(
    IWS_mds %>% mutate(dataset = "IWS", mol_manual = MOLECULAR_GROUP) %>% select(IPSSM, mol_manual, dataset), 
    gesmd_dataset %>% mutate(dataset = "GESMD") %>% select(IPSSM, mol_manual, dataset)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
    Taxonomy = case_when(
        mol_manual == "\"-7/SETBP1\"" ~ "7-/SETBP1",
        mol_manual %in% c("mNOS", "No-event") ~ "Other",
        TRUE ~ mol_manual
    )) %>%
  filter(!is.na(IPSSM)) %>%
  filter(Taxonomy %in% c("bi-TET2", "EZH2-ASXL1", "7-/SETBP1", "IDH-STAG2", "Other")) %>%
  mutate(Taxonomy = factor(Taxonomy, levels = c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other"))) %>%
  group_by(Taxonomy, IPSSM, dataset) %>%
  summarize(N = n()) %>%
  group_by(Taxonomy, dataset) %>%
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
    x = "Prognostic group",
    fill = "IPSSM") +
    facet_grid(~ Taxonomy) +
  theme(plot.title = element_text(hjust = 0.5)) 
png("figures/GESMD_IWS_clustering/taxonomy_comparison/IPSSM_summary.png", res = 300, height = 1200, width = 1800)
ipssm_summary
dev.off()


## Prognosis
joint_prognosis_plot <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS", mol_manual = MOLECULAR_GROUP) %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), IPSSM, IPSSM_SCORE, mol_manual, sub_group, AGE, SEX, dataset, BM_BLAST),
    gesmd_dataset %>% mutate(dataset = "GESMD") %>% 
    select(ends_with("STATUS"), ends_with("YEARS"), IPSSM, IPSSM_SCORE,  mol_manual, sub_group, AGE, SEX, dataset, BM_BLAST)
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
    Taxonomy = case_when(
        mol_manual == "\"-7/SETBP1\"" ~ "7-/SETBP1",
        mol_manual %in% c("mNOS", "No-event") ~ "Other",
        TRUE ~ mol_manual
    ),
    Taxonomy = factor(Taxonomy),
    IPSSM = factor(IPSSM, levels = c( "Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High"))) 


joint_prognosis <- joint_prognosis_plot %>% 
  mutate(Taxonomy = relevel(Taxonomy, ref = "Other"),
  IPSSM = factor(IPSSM, levels = c("Low", "Very-Low", "Moderate-Low", "Moderate-High", "High", "Very-High")))


taxons <- c("bi-TET2", "EZH2-ASXL1", "7-/SETBP1", "IDH-STAG2", "Other")
main_model_tabs <- lapply(taxons, function(group){

  joint_prognosis_sub <- joint_prognosis %>%
    filter(Taxonomy == group) %>%
    mutate(Taxonomy = droplevels(Taxonomy),
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

colors6 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#999999", "grey40",  "black")

hr_os_group_plot <- bind_rows(main_model_tabs,
  tibble(Variable = "Moderate-High", Group = taxons, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other"))) %>%
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

png("figures/GESMD_IWS_clustering/taxonomy_comparison/OS_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_os_group_plot
dev.off()


joint_prognosis_os_plot <- joint_prognosis %>%
    mutate(comb_group = paste(Taxonomy, IPSSM, sep = "_"),
    comb_group = factor(comb_group, levels = unique(comb_group)),
    comb_group = relevel(comb_group, ref = "Other_Low"))

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
    Variable = "Other_Low",
    Group = "Other",
    IPSSM = "Low",
    HR = 1,
    HR_Low = 1,
    HR_High = 1,
    p_value = NA
  ))



hr_os_group_plot2 <- os_mod_subgroup_ipssm_tab %>%
    filter(Group %in% c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other")) %>%
mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other"))) %>%
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
  ggtitle("OS (Ref: Other/IPSSM Low)")

png("figures/GESMD_IWS_clustering/taxonomy_comparison/OS_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
hr_os_group_plot2
dev.off()

png("figures/GESMD_IWS_clustering/taxonomy_comparison/OS_subgroups_HRs_panel.png", width = 2500, height = 2000, res = 300)
plot_grid(hr_os_group_plot2, hr_os_group_plot, ncol = 1, labels = c("A", "B"))
dev.off()



## AMLt
IWS_aml <- filter(joint_prognosis, dataset == "IWS")

aml_model_tabs <- lapply(taxons, function(group){

  IWS_sub <- IWS_aml %>%
    filter(Taxonomy == group) %>%
    mutate(Taxonomy = droplevels(Taxonomy),
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
  tibble(Variable = "Moderate-High", Group = taxons, HR = 1, HR_Low = 1, HR_High = 1, p_value = NA)
)  %>%
mutate(Variable = factor(Variable, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other"))) %>%
      filter(HR_Low != 0) %>%
     filter(Group != "EZH2-ASXL1") %>%
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

png("figures/GESMD_IWS_clustering/taxonomy_comparison/AMLt_subgroups_HRs.png", width = 2500, height = 1200, res = 300)
hr_aml_group_plot
dev.off()


IWS_aml_plot <- IWS_aml %>%
    mutate(comb_group = paste(Taxonomy, IPSSM, sep = "_"),
    comb_group = factor(comb_group, levels = unique(comb_group)),
    comb_group = relevel(comb_group, ref = "Other_Low"))

cat_count_aml <- table(IWS_aml_plot$comb_group)
sel_cats_aml <- names(cat_count_aml[cat_count_aml >= 10])
sel_cats_aml <- sel_cats_aml[grepl("EZH2|SETBP1|bi-TET2|STAG2|Other", sel_cats_aml)]
sel_cats_aml <- sel_cats_aml[!sel_cats_aml %in% c("EZH2-ASXL1_Moderate-High", "Other_NA", "bi-TET2_Very-Low", "Other_Very-Low")]
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
    Variable = "Other_Low",
    Group = "Other",
    IPSSM = "Low",
    HR = 1,
    HR_Low = 1,
    HR_High = 1,
    p_value = NA
  ))


hr_aml_group_plot2 <- aml_mod_subgroup_ipssm_tab %>%
mutate(IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
      Group = factor(Group, levels = c("EZH2-ASXL1", "bi-TET2", "7-/SETBP1", "IDH-STAG2", "Other"))) %>%
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
  ggtitle("AMLt (Ref: Other/IPSSM Low)")

png("figures/GESMD_IWS_clustering/taxonomy_comparison/AMLt_subgroups_HRs_relative.png", width = 2500, height = 1200, res = 300)
hr_aml_group_plot2
dev.off()


png("figures/GESMD_IWS_clustering/taxonomy_comparison/AMLt_subgroups_HRs_panel.png", width = 2500, height = 2000, res = 300)
plot_grid(hr_aml_group_plot2, hr_aml_group_plot, ncol = 1, labels = c("A", "B"))
dev.off()


main_os_mol <- coxph(Surv(OS_YEARS,OS_STATUS) ~ mol_manual + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) 
int_os_mol <- coxph(Surv(OS_YEARS,OS_STATUS) ~ mol_manual*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) 


coxph(Surv(OS_YEARS,OS_STATUS) ~ mol_manual + IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) 


coxph(Surv(OS_YEARS,OS_STATUS) ~ mol_manual*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) 
coxph(Surv(OS_YEARS,OS_STATUS) ~ sub_group*IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis) 


coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = mol_manual == "IDH/STAG2") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = joint_prognosis, subset = sub_group == "STAG2") 


coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = mol_manual == "IDH/STAG2" & dataset == "IWS") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = mol_manual == "IDH/STAG2"  & dataset == "GESMD") 

coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = sub_group == "STAG2" & dataset == "IWS") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = joint_prognosis, subset = sub_group == "STAG2"  & dataset == "GESMD") 

coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX, data = joint_prognosis) 


jp_stag2 <- joint_prognosis %>%
    filter(IPSSM %in% c("Moderate-High", "High", "Very-High")) %>%
    mutate(IPSSM = droplevels(IPSSM))

coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX + dataset, data = jp_stag2, subset = mol_manual == "IDH/STAG2") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX + dataset, data = jp_stag2, subset = sub_group == "STAG2") 

coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = jp_stag2, subset = mol_manual == "IDH/STAG2" & dataset == "IWS") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM + AGE + SEX, data = jp_stag2, subset = sub_group == "STAG2" & dataset == "IWS") 

coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = jp_stag2, subset = mol_manual == "IDH/STAG2") 
coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSM_SCORE + AGE + SEX + dataset, data = jp_stag2, subset = sub_group == "STAG2") 



