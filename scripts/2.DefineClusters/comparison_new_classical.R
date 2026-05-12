#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compare classical subgroups with current versions
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for subgroups exploration
#' 
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.9 R
#'
#' ---------------------------

# Load libraries and data
library(MASS)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(survminer)
library(survival)

load("results/GESMD_IWS_clustering/gesmd_IWS_mds.Rdata")
load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")
load("results/hershberger/hershberger_mds.Rdata")
load("results/hershberger/hershberger_full.Rdata")



joint_full <- bind_rows(
    IWS_full %>% mutate(dataset = "IWS") %>% 
    mutate(complex = ifelse(complex == "complex", 1, 0),
    CYTO_IPSSR = case_when(
        CYTO_IPSSR == "Very-Poor" ~ "Very Poor",
        CYTO_IPSSR == "Very-Good" ~ "Very Good",
        TRUE ~ CYTO_IPSSR
    )
    ),
    gesmd_full %>% mutate(dataset = "GESMD") %>%
    mutate(CYTO_IPSSR = case_when(
        CYTO_IPSSR == "Intermediate" ~ "Int",
        TRUE ~ CYTO_IPSSR
    )
    ),
    hersh_all_mds %>% mutate(dataset = "MLL") %>%
    mutate(CYTO_IPSSR = case_when(
        CYTO_IPSSR == "Intermediate" ~ "Int",
        TRUE ~ CYTO_IPSSR
    ),
    consensus = WHO
)
) %>% mutate(CYTO_IPSSR = factor(CYTO_IPSSR, levels = c( "Very Good", "Good",  "Int", "Poor", "Very Poor")))


joint_mds <- bind_rows(
    IWS_mds %>% mutate(dataset = "IWS"),
    gesmd_dataset %>% mutate(dataset = "GESMD"),
    hersh_mds %>% mutate(dataset = "MLL") 
) 

joint_full_subgroup <- left_join(select(joint_full, -sub_group) , select(joint_mds, ID, sub_group), by = "ID") %>%
    mutate(sub_group_comb = ifelse(is.na(sub_group), consensus, as.character(sub_group)),
    sub_group_comb = case_when(
        sub_group_comb %in% c("mutated SF3B1", "MDS-SF3B1") ~ "SF3B1",
        sub_group_comb == "MDS-5q" ~ "del5q",
        sub_group_comb %in% c("Mutated TP53", "MDS-TP53") ~ "TP53",
        sub_group_comb == "Low blasts" ~ "MDS-LB",
        TRUE ~ sub_group_comb
    ))

## Define dataset for comparison
groups <- c("TP53", "Complex", "del5q", "del5q-IB", "SF3B1", "SF3B1-IB",  "MDS-IB1", "MDS-IB2")

joint_full_filter <- filter(joint_full_subgroup, sub_group_comb %in% groups) 
joint_full_filter$sub_group_comb <- factor(joint_full_filter$sub_group_comb, levels = groups)

clin_joint <- filter(joint_full_filter, dataset %in% c("IWS", "GESMD"))
hersh_filt <- filter(joint_full_filter, dataset == "MLL") %>%
    select(-ANC, -MONOCYTES) 

colors <- c("#FF9E59", "#7A3500", "#FFF9AE", "#B8A600", "#56B4E9", "#003E61",  "grey40",  "black")

## Clinical variables
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")



## Test for negative binomial (BM_BLAST)
res_bm <- lapply(seq_len(length(groups)), function(i){
      cl <- groups[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group_comb, cl))
     poiss_lm <- summary(glm.nb(BM_BLAST ~ cluster + AGE + SEX + dataset, tmp))
     coefs <- poiss_lm$coefficients[-1, ]

     tmp_mll <- mutate(hersh_filt, cluster = relevel(sub_group_comb, cl))
     poiss_lm_mll <- summary(glm.nb(BM_BLAST ~ cluster + AGE + SEX, tmp_mll))
    coefs_mll <- poiss_lm_mll$coefficients[-1, ]

     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
              filter(Comp_clust %in% groups) %>%
              left_join(as_tibble(coefs_mll[, c(1, 4)]) %>%
                          mutate(Ref_cluster = cl, 
                                 Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                                 Stat = exp(Estimate)) %>%
                          filter(Comp_clust %in% groups),
            by = c("Ref_cluster", "Comp_clust"), suffix = c("", "_MLL")) 

     colnames(coef_df)[c(2, 7)] <- c("P_value", "P_value_MLL")
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value, Stat_MLL, P_value_MLL)
    }) %>%
    Reduce(., f = rbind) %>%
     mutate(Clin_var = "BM_BLAST") %>%
    filter(Stat > 1) %>%
    arrange(P_value)
  

poisson <- c("WBC", "ANC", "MONOCYTES")
res_poisson <- lapply(poisson, function(var){
    cell_test <- lapply(seq_len(length(groups) ), function(i){
      cl <- groups[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group_comb, cl))
     poiss_lm <- summary(glm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp, 
                      family = "poisson"))
     coefs <- poiss_lm$coefficients[-1, ]

    if (var %in% colnames(hersh_filt)) {

     tmp_mll <- mutate(hersh_filt, cluster = relevel(sub_group_comb, cl))
     poiss_lm_mll <- summary(glm(formula (paste(var, " ~ cluster + AGE + SEX ")), tmp_mll,
             family = "poisson"))
     coefs_mll <- poiss_lm_mll$coefficients[-1, ]
     mll_tib <- as_tibble(coefs_mll[, c(1, 4)])
    } else {
       coefs_mll <-  mll_tib <- data.frame(Estimate = NA, `Pr(>|z|)` = NA, cluster = groups[groups != cl])
       rownames(coefs_mll) <- paste0("cluster", groups[groups != cl])
        mll_tib <- select(mll_tib, -cluster)
    }
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
        filter(Comp_clust %in% groups) %>%
        left_join(mll_tib %>%
                 mutate(Ref_cluster = cl,
                     Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                     Stat = exp(Estimate)) %>%
                 filter(Comp_clust %in% groups),
         by = c("Ref_cluster", "Comp_clust"), suffix = c("", "_MLL"))

     colnames(coef_df)[c(2, 7)] <- c("P_value", "P_value_MLL")
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value, Stat_MLL, P_value_MLL)
    }) %>%
    Reduce(., f = rbind) 
  out <- mutate(cell_test, 
                Clin_var = var) %>%
    filter(Stat > 1) %>%
    arrange(P_value)
  out
})
lapply(res_poisson, function(x) filter(x, Ref_cluster %in% groups & Comp_clust %in% groups))



normal <- c("HB", "PLT")
res_norm <- lapply(normal, function(var){
    cell_test <- lapply(seq_len(length(groups)), function(i){
      cl <- groups[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group_comb, cl))
     lm <- summary(lm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp))
     coefs <- lm$coefficients[-1, ]


     tmp_mll <- mutate(hersh_filt, cluster = relevel(sub_group_comb, cl))
     lm_mll <- summary(lm(formula (paste(var, " ~ cluster + AGE + SEX ")), tmp_mll))
     coefs_mll <- lm_mll$coefficients[-1, ]
     mll_tib <- as_tibble(coefs_mll[, c(1, 4)])


     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = Estimate) %>%
        filter(Comp_clust %in% groups) %>%
        left_join(mll_tib %>%
                 mutate(Ref_cluster = cl,
                     Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                     Stat = Estimate) %>%
                 filter(Comp_clust %in% groups),
         by = c("Ref_cluster", "Comp_clust"), suffix = c("", "_MLL"))
     colnames(coef_df)[c(2, 7)] <- c("P_value", "P_value_MLL")
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value, Stat_MLL, P_value_MLL)
    }) %>%
    Reduce(., f = rbind) 
  out <- mutate(cell_test, 
                Clin_var = var) %>%
    filter(Stat > 0) %>%
    arrange(P_value)
  out
})

clin_test <- rbind(res_bm,
                   Reduce(res_poisson, f = rbind),
                   Reduce(res_norm, f = rbind)) 

write.table(clin_test, 
            file = "results/GESMD_IWS_clustering/classical_subgroup_clinical_tests.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)


colores_blasts <- c("#FFECB3", "#FFB300", "#D32F2F")
colores_wbc <- c("#E1BEE7", "#BA68C8", "#8E24AA", "#4A148C")
colores_monocitos <- c("#CFD8DC", "#4FC3F7", "#FB8C00")
colores_hb <- c("#F8BBD0", "#F06292", "#E91E63", "#880E4F")
colores_plt <- c("#C8E6C9", "#81C784", "#43A047", "#1B5E20")


makeBarPlot <- function(df, var){
    df$dataset <- fct_rev(df$dataset)
   ggplot(df, aes(x = dataset, fill = Cell_class)) +
    geom_bar(position = "fill") +
    facet_grid(sub_group_comb ~ .) +
    coord_flip() +
    theme_bw() +
    guides(fill = guide_legend(ncol = 2)) +
    theme(legend.position = "top",
    legend.direction = "vertical",
    axis.text.y = element_blank(),
         strip.text = element_blank(),
          strip.background = element_blank(),
        plot.margin = margin(r = 1, l = 1, unit = "pt")) +
    xlab("") +
    scale_y_continuous(name = "", labels = scales::label_percent())
}

## BM blasts
blast_plot <- joint_full_filter %>%
    mutate(Cell_class = case_when(
        BM_BLAST <= 5 ~ "<5%",
        BM_BLAST <= 10 ~ "5-10%",
        TRUE ~ ">=10%"
    ),
    Cell_class = factor(Cell_class, levels = c("<5%", "5-10%", ">=10%"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "BM Blasts", values = colores_blasts)

## WBC
wbc_plot <- joint_full_filter %>%
    mutate(Cell_class = case_when(
        WBC < 1.5 ~ "<1.5",
        WBC < 3 ~ "1.5-3",
        WBC < 4 ~ "3-4",
        TRUE ~ ">=4"
    ),
    Cell_class = factor(Cell_class, levels = c("<1.5", "1.5-3", "3-4", ">=4"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "WBC", values = colores_wbc)

## MONO
mono_plot <- joint_full_filter %>%
    filter(dataset != "MLL") %>%
    mutate(Cell_class = case_when(
        MONOCYTES < 0.5 ~ "<0.5",
        MONOCYTES < 1 ~ "0.5-1",
        TRUE ~ ">=1"
    ),
    Cell_class = factor(Cell_class, levels = c("<0.5", "0.5-1", ">=1"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "MONOCYTES", values = colores_monocitos)

## ANC
anc_plot <- joint_full_filter %>%
    filter(dataset != "MLL") %>%    
    mutate(Cell_class = case_when(
        ANC < 0.4 ~ "<0.4",
        ANC < 1 ~ "0.4-1",
        ANC < 1.5 ~ "1-1.5",
        TRUE ~ ">=1.5"
    ),
    Cell_class = factor(Cell_class, levels = c("<0.4", "0.4-1", "1-1.5", ">=1.5"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "ANC", values = colores_wbc)


## HB
hb_plot <- joint_full_filter %>%
    mutate(Cell_class = case_when(
        HB < 8 ~ "<8",
        HB < 10 ~ "8-10",
        HB < 12 ~ "10-12",
        TRUE ~ ">=12"
    ),
    Cell_class = factor(Cell_class, levels = c("<8", "8-10", "10-12", ">=12"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "HB", values = colores_hb)

## PLT
plt_plot <- joint_full_filter %>%
    mutate(Cell_class = case_when(
        PLT < 50 ~ "<50",
        PLT < 100 ~ "50-100",
        PLT < 150 ~ "100-150",
        TRUE ~ ">=150"
    ),
    Cell_class = factor(Cell_class, levels = c("<50", "50-100", "100-150", ">=150"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "PLT", values = colores_plt)

clin_panel <- plot_grid(blast_plot + theme(axis.text.y = element_text()), 
    wbc_plot , 
    hb_plot , 
    plt_plot + theme(strip.text = element_text(), strip.background = element_rect()), 
    nrow = 1)
png("figures/GESMD_IWS_clustering/classical_groups_exploration/clin_vars_subgroups_barplot.png", width = 3000, height = 2300, res = 300)
clin_panel
dev.off()

clin_panel_sup <- plot_grid(mono_plot + theme(axis.text.y = element_text()), 
    anc_plot + theme(strip.text = element_text(), strip.background = element_rect()), 
    nrow = 1)

png("figures/GESMD_IWS_clustering/classical_groups_exploration/clin_vars_subgroups_barplot_sup.png", width = 2000, height = 2300, res = 300)
clin_panel_sup
dev.off()


## Molecular
# Select mutations Freq > 4% and not included in sub-groups
mutations <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "PHF6", "plus8", "delY", "del20q", "del7")

hersh_filt$TET2other <- hersh_filt$TET2mono
mut_ORs <- lapply(mutations, function(variable){
    lapply(seq_len(length(groups) ), function(i){
        ref_cl <- groups[i]
        lapply(groups[groups != ref_cl], function(comp_cl){
            
            sel_df <- filter(clin_joint, sub_group_comb %in% c(ref_cl, comp_cl)) %>%
                mutate(cluster = factor(sub_group_comb, levels = c(ref_cl, comp_cl)))
            tab <- table(sel_df[[variable]], sel_df$cluster)
            fisher_res <- tryCatch(fisher.test(tab), error = function(e) NULL)

            sel_df_mll <- filter(hersh_filt, sub_group_comb %in% c(ref_cl, comp_cl)) %>%
                mutate(cluster = factor(sub_group_comb, levels = c(ref_cl, comp_cl)))
            tab_mll <- table(sel_df_mll[[variable]], sel_df_mll$cluster)
            fisher_res_mll <- tryCatch(fisher.test(tab_mll), error = function(e) NULL)

            tibble(Gene = variable,
                   Ref_cluster = ref_cl,
                   Comp_clust = comp_cl,
                   OR = ifelse(is.null(fisher_res), NA, fisher_res$estimate),
                   P_value = ifelse(is.null(fisher_res), NA, fisher_res$p.value),
                   OR_MLL = ifelse(is.null(fisher_res_mll), NA, fisher_res_mll$estimate),
                   P_value_MLL = ifelse(is.null(fisher_res_mll), NA, fisher_res_mll$p.value))
        }) %>% Reduce(., f = rbind)
    }) %>% Reduce(., f = rbind)
}) %>% Reduce(., f = rbind) %>%
    filter(OR > 1) %>%
    arrange(Gene, P_value)


write.table(mut_ORs, 
            file = "results/GESMD_IWS_clustering/subgroup_mutation_tests_classical.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

mut_gene_plot <- c(
  "ASXL1", "SRSF2", "DNMT3A", "TET2other", 
  "RUNX1", "U2AF1", "BCOR", "ZRSR2", 
  "IDH2", "SETBP1", "DDX41", "PHF6",
  "plus8", "delY", "del20q", "del7"
)

mut_tib <- joint_full_filter %>%
    mutate(TET2other = ifelse(is.na(TET2other), TET2mono, TET2other),
        sub_group_comb = droplevels(sub_group_comb)) %>%
    group_by(sub_group_comb, dataset) %>%
    summarize_at(vars(all_of(mut_gene_plot)), mean, na.rm = TRUE) 

mut_mat <- as.matrix(mut_tib[, -c(1,2)])
rownames(mut_mat) <- seq_len(nrow(mut_tib))
colnames(mut_mat)[colnames(mut_mat) == "TET2other"] <- "TET2-mono"
colnames(mut_mat)[colnames(mut_mat) == "plus8"] <- "+8"
colnames(mut_mat)[colnames(mut_mat) == "delY"] <- "-Y"
colnames(mut_mat)[colnames(mut_mat) == "del7"] <- "-7"

row_annot <- data.frame(mut_tib[, 1])

mut_heatmap <- Heatmap(mut_mat, 
        row_labels = mut_tib$dataset,
        show_heatmap_legend = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,     
        row_names_side = "left",
        column_names_side = "top",
        cluster_rows = FALSE, 
        cluster_columns = TRUE, 
        show_column_dend = FALSE, 
        row_split = row_annot$sub_group, 
        row_title = paste0("\n", unique(row_annot$sub_group)), 
        row_gap = unit(2, "mm"), 
        row_title_gp = gpar(
          col = c("black", "white", "black", "black", "black", "white", "black", "white"), 
          fill = colors
        ),        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", mut_mat[i, j] * 100), x, y, 
                    gp = gpar(fontsize = 10))
        },
        border = TRUE
)

png("figures/GESMD_IWS_clustering/classical_groups_exploration/mutation_subgroups_pheatmap.png", width = 2000, height = 2300, res = 300)
mut_heatmap
dev.off()


ipssm_summary <- joint_full_filter %>% filter(!is.na(IPSSM) ) %>%
  group_by(sub_group_comb, IPSSM, dataset) %>%
  summarize(N = n()) %>%
  group_by(sub_group_comb, dataset) %>%
  mutate(Freq = N/sum(N),
         IPSSM = fct_recode(IPSSM, 
         "VL" = "Very-Low",
         "VL" = "Very Low",
         "L" = "Low",
         "ML" = "Moderate-Low", 
         "ML" = "Moderate Low", 
         "MH" = "Moderate-High",
         "MH" = "Moderate High",
         "H" = "High",
         "VH" = "Very-High",
         "VH" = "Very High"),
         IPSSM = factor(IPSSM, levels = c("VL", "L", "ML", "MH", "H", "VH"))) %>%
  ggplot(aes(x = fct_rev(dataset), y = Freq*100, fill = fct_rev(IPSSM))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "IPSSM classification",
    y = "",
    x = "",
    fill = "IPSSM") +
    facet_grid(sub_group_comb ~ .) +
   coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
  legend.position = "top", strip.background = element_rect()) 
png("figures/GESMD_IWS_clustering/classical_groups_exploration/IPSSM_summary.png", res = 300, height = 900, width = 3000)
ipssm_summary
dev.off()



age_plot <-  ggplot(joint_full_filter, aes(x = dataset, y = AGE, fill = sub_group_comb)) +
      geom_boxplot() +
      facet_grid(~ sub_group_comb) +
         scale_fill_manual(values = colors, name = "Sub-group") +
      ylab("Age (Years)") +
      xlab("") +
      ggtitle("Age") +
      coord_cartesian(ylim = c(35, 100)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))

png("figures/GESMD_IWS_clustering/classical_groups_exploration/age.png", res = 300, height = 1400, width = 2000)
age_plot
dev.off()

summary(lm(AGE ~ sub_group_comb + dataset, mutate(clin_joint, sub_group_comb = relevel(sub_group_comb, "TP53"))))
summary(lm(AGE ~ sub_group_comb + dataset, mutate(clin_joint, sub_group_comb = relevel(sub_group_comb, "del5q"))))
summary(lm(AGE ~ sub_group_comb + dataset, mutate(clin_joint, sub_group_comb = relevel(sub_group_comb, "SF3B1"))))
summary(lm(AGE ~ sub_group_comb + dataset, mutate(clin_joint, sub_group_comb = relevel(sub_group_comb, "MDS-IB1"))))



ht_grob <- grid.grabExpr(draw(mut_heatmap))
png("figures/GESMD_IWS_clustering/classical_groups_exploration/panel_summary.png", res = 300, height = 4700, width = 3500)
plot_grid(
    plot_grid(clin_panel, ipssm_summary, ncol = 2, rel_widths = c(2, 1), labels = c("A", "C")),
    ht_grob, ncol = 1, labels = c("", "B"), rel_heights = c(1, 1)
)
dev.off()




surv_all <- survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sub_group_comb, joint_full_filter) %>%
    ggsurvplot(data = joint_full_filter, surv.median.line = "hv", palette = colors,
     risk.table = TRUE, break.time.by = 2, xlim = c(0, 10),
      legend.labs  = levels(joint_full_filter$sub_group_comb))
png("figures/GESMD_IWS_clustering/classical_groups_exploration/OS_subgroups_survival.png", width = 1500, height = 2000, res = 300)
surv_all$plot + theme(plot.title = element_text(hjust = 0.5)) + labs(title = "Overall Survival by sub-group", x = "Time (years)", y = "Survival probability", color = "Sub-group")
dev.off()


joint_full_filter %>%
    filter(sub_group_comb %in% c("SF3B1", "SF3B1-IB", "MDS-IB1", "MDS-IB2")) %>%
    mutate(cl = droplevels(sub_group_comb)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ cl*IPSSM_SCORE + AGE + SEX, .)

joint_full_filter %>%
    filter(sub_group_comb %in% c("del5q", "del5q-IB", "MDS-IB1", "MDS-IB2")) %>%
    mutate(cl = droplevels(sub_group_comb)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ cl* IPSSM_SCORE + AGE + SEX, .)

joint_full_filter %>%
    filter(sub_group_comb %in% c("TP53", "Complex", "MDS-IB1", "MDS-IB2")) %>%
    mutate(cl = droplevels(sub_group_comb)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ cl* IPSSM_SCORE + AGE + SEX, .)



## Plot median OS
ipssm_raw_surv <- survfit( Surv(OS_YEARS,OS_STATUS) ~ sub_group_comb+IPSSM, clin_joint)
ipssm_raw_surv_tib <- data.frame(summary(ipssm_raw_surv)$table) %>%
  rownames_to_column("Group") %>%
  separate(Group, into = c("sub_group", "IPSSM"), sep = "\\, ") %>%
  mutate(sub_group = gsub("sub_group_comb=", "", sub_group),
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

png("figures/GESMD_IWS_clustering/classical_groups_exploration/OS_subgroups_median.png", width = 1500, height = 2000, res = 300)
median_OS_plot
dev.off()


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
    df <- subset(clin_joint, sub_group_comb == cl)
    df_est <- getEstimates(gene, df, joint = TRUE) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind)
}) %>% Reduce(f = rbind)
mut_hr_full_joint <- lapply(test_muts, function(gene){
  df_est <- getEstimates(gene, clin_joint, joint = TRUE) %>%
      mutate(sub_group = "Both cohorts")
  }) %>% Reduce(f = rbind)


cyto_hr <- lapply(groups, function(cl){
    df <- subset(clin_joint, sub_group_comb == cl)
    df_est <- getEstimates("as.numeric(CYTO_IPSSR)", df, joint = TRUE) %>%
      mutate(sub_group = cl)
  }) %>% Reduce(f = rbind) %>%
  bind_rows(., getEstimates("as.numeric(CYTO_IPSSR)", clin_joint, joint = TRUE) %>%
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


