#' ---------------------------
#'
#' Purpose of script:
#'
#'  Explore the subgroups of MDS
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
library(ggh4x)
library(cowplot)
library(ComplexHeatmap)

load("results/GESMD_IWS_clustering/gesmd_IWS_mds.Rdata")
load("results/hershberger/hershberger_mds.Rdata")


# colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
#     "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")
colors_all <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
    "#D55E00", "#999999", "grey40",  "black")
colors <-  c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
     "#999999", "grey40",  "black")
scale_col <- scale_color_manual(values = colors, name = "Sub-group")
scale_fill <- scale_fill_manual(values = colors, name = "Sub-group")

## Proportion of each subgroup in each dataset
mean(IWS_mds$sub_group %in% c("EZH2", "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex"))
mean(gesmd_dataset$sub_group %in% c("EZH2", "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex"))
mean(hersh_mds$sub_group %in% c("EZH2", "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex"))

prop.table(table(IWS_mds$sub_group))
prop.table(table(gesmd_dataset$sub_group))
prop.table(table(hersh_mds$sub_group))

mean(IWS_mds$sub_group %in% c("del5q-IB", "SF3B1-IB", "Complex"))
mean(gesmd_dataset$sub_group %in% c("del5q-IB", "SF3B1-IB", "Complex"))
mean(hersh_mds$sub_group %in% c("del5q-IB", "SF3B1-IB", "Complex"))


## Clinical variables
clusters <- levels(IWS_mds$sub_group)
clusters <- c("EZH2", "TET2-bi", "-7", "STAG2", "MDS-LB", "MDS-IB1", "MDS-IB2")



## Tests
sub_groups <- c("EZH2", "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex", "MDS-LB", "MDS-IB1", "MDS-IB2")
sub_groups <- c("EZH2", "TET2-bi", "-7", "STAG2", "MDS-LB", "MDS-IB1", "MDS-IB2")

clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_joint3 <- bind_rows(
    IWS_mds %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "GESMD"),
    hersh_mds %>% select(any_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "MLL")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD", "MLL"))) %>%
    filter(sub_group %in% sub_groups)

clin_joint <- clin_joint3 %>%
    filter(dataset != "MLL") %>%
    mutate(sub_group = droplevels(sub_group),
    dataset = droplevels(dataset))

## Compute cytopenia prop
#nb <- c("BM_BLAST", "PB_BLAST")
nb <- c("BM_BLAST") ## Exclude PB BLAST

res_nb <- lapply(nb, function(var){
    cell_test <- lapply(seq_len(length(clusters)), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm.nb(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp))
     coefs <- poiss_lm$coefficients[-1, ]

     tmp_mll <- mutate(hersh_mds, cluster = relevel(sub_group, cl))
     poiss_lm_mll <- summary(glm.nb(formula (paste(var, " ~ cluster + AGE + SEX ")), tmp_mll))
    coefs_mll <- poiss_lm_mll$coefficients[-1, ]

     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
              filter(Comp_clust %in% sub_groups) %>%
              left_join(as_tibble(coefs_mll[, c(1, 4)]) %>%
                          mutate(Ref_cluster = cl, 
                                 Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                                 Stat = exp(Estimate)) %>%
                          filter(Comp_clust %in% sub_groups),
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
lapply(res_nb, function(x) filter(x, Ref_cluster %in% sub_groups & Comp_clust %in% sub_groups))


poisson <- c("WBC", "ANC", "MONOCYTES")
res_poisson <- lapply(poisson, function(var){
    cell_test <- lapply(seq_len(length(clusters) ), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp, 
                      family = "poisson"))
     coefs <- poiss_lm$coefficients[-1, ]

    if (var %in% colnames(hersh_mds)) {

     tmp_mll <- mutate(hersh_mds, cluster = relevel(sub_group, cl))
     poiss_lm_mll <- summary(glm(formula (paste(var, " ~ cluster + AGE + SEX ")), tmp_mll,
             family = "poisson"))
     coefs_mll <- poiss_lm_mll$coefficients[-1, ]
     mll_tib <- as_tibble(coefs_mll[, c(1, 4)])
    } else {
       coefs_mll <-  mll_tib <- data.frame(Estimate = NA, `Pr(>|z|)` = NA, cluster = sub_groups[sub_groups != cl])
       rownames(coefs_mll) <- paste0("cluster", sub_groups[sub_groups != cl])
        mll_tib <- select(mll_tib, -cluster)
    }
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
        filter(Comp_clust %in% sub_groups) %>%
        left_join(mll_tib %>%
                 mutate(Ref_cluster = cl,
                     Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                     Stat = exp(Estimate)) %>%
                 filter(Comp_clust %in% sub_groups),
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
lapply(res_poisson, function(x) filter(x, Ref_cluster %in% sub_groups & Comp_clust %in% sub_groups))

normal <- c("HB", "PLT")
res_norm <- lapply(normal, function(var){
    cell_test <- lapply(seq_len(length(clusters)), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     lm <- summary(lm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp))
     coefs <- lm$coefficients[-1, ]


     tmp_mll <- mutate(hersh_mds, cluster = relevel(sub_group, cl))
     lm_mll <- summary(lm(formula (paste(var, " ~ cluster + AGE + SEX ")), tmp_mll))
     coefs_mll <- lm_mll$coefficients[-1, ]
     mll_tib <- as_tibble(coefs_mll[, c(1, 4)])


     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = Estimate) %>%
        filter(Comp_clust %in% sub_groups) %>%
        left_join(mll_tib %>%
                 mutate(Ref_cluster = cl,
                     Comp_clust = gsub("cluster", "", rownames(coefs_mll)),
                     Stat = Estimate) %>%
                 filter(Comp_clust %in% sub_groups),
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

clin_test <- rbind(Reduce(res_nb, f = rbind),
                   Reduce(res_poisson, f = rbind),
                   Reduce(res_norm, f = rbind)) 
    
write.table(clin_test, 
            file = "results/GESMD_IWS_clustering/subgroup_clinical_tests.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
subset(clin_test, Ref_cluster %in% sub_groups & Comp_clust %in% sub_groups & P_value < 0.05)


colores_blasts <- c("#FFECB3", "#FFB300", "#D32F2F")
colores_wbc <- c("#E1BEE7", "#BA68C8", "#8E24AA", "#4A148C")
colores_monocitos <- c("#CFD8DC", "#4FC3F7", "#FB8C00")
colores_hb <- c("#F8BBD0", "#F06292", "#E91E63", "#880E4F")
colores_plt <- c("#C8E6C9", "#81C784", "#43A047", "#1B5E20")


makeBarPlot <- function(df, var){
    df$dataset <- fct_rev(df$dataset)

    df2 <- group_by(df, sub_group, dataset) %>% 
        mutate(dataset2 = paste0(dataset, " (", n(), ")"))

    levs <- arrange(df2, desc(dataset)) %>% pull(dataset2) %>% unique()
    df2$dataset2 <- factor(df2$dataset2, levels = levs)

    ggplot(df2, aes(y = dataset2, fill = Cell_class)) +
        geom_bar(position = "fill", orientation = "y") + 
        facet_grid2(sub_group ~ ., 
                    switch = "y", 
                    scales = "free_y",
                    independent = "y",
                    strip = strip_themed(
                        background_y = elem_list_rect(fill = colors)
                ))+
        theme_bw() +
        guides(fill = guide_legend(ncol = 2)) +
        theme(
            legend.position = "top",
            legend.direction = "vertical",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.y = element_blank(),
            plot.margin = margin(r = 1, l = 1, unit = "pt")
        ) +
        xlab("") +
        scale_x_continuous(name = "", labels = scales::label_percent()) +
        scale_y_discrete(name = "", position = "right")
}

## BM blasts
blast_plot <- clin_joint3 %>%
    mutate(Cell_class = case_when(
        BM_BLAST <= 5 ~ "<5%",
        BM_BLAST <= 10 ~ "5-10%",
        TRUE ~ ">=10%"
    ),
    Cell_class = factor(Cell_class, levels = c("<5%", "5-10%", ">=10%"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "BM Blasts", values = colores_blasts)

## WBC
wbc_plot <- clin_joint3 %>%
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
mono_plot <- clin_joint3 %>%
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
anc_plot <- clin_joint3 %>%
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
hb_plot <- clin_joint3 %>%
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
plt_plot <- clin_joint3 %>%
    mutate(Cell_class = case_when(
        PLT < 50 ~ "<50",
        PLT < 100 ~ "50-100",
        PLT < 150 ~ "100-150",
        TRUE ~ ">=150"
    ),
    Cell_class = factor(Cell_class, levels = c("<50", "50-100", "100-150", ">=150"))) %>%
    makeBarPlot() +
    scale_fill_manual(name = "PLT", values = colores_plt)
 
clin_panel <- plot_grid(blast_plot + theme(strip.text = element_text(),
    strip.text.y.left = element_text(angle = 0, color = "white", face = "bold")),
    wbc_plot , 
    hb_plot , 
    plt_plot + theme(axis.text.y = element_text(), axis.ticks.y = element_line()), 
    nrow = 1, rel_widths = c(1.2, 1, 1, 1.5))
png("figures/GESMD_IWS_clustering/subgroup_exploration/clin_vars_subgroups_barplot.png", width = 3000, height = 2000, res = 300)
clin_panel
dev.off()

clin_panel_sup <- plot_grid(mono_plot + theme(strip.text = element_text(),
    strip.text.y.left = element_text(angle = 0, color = "white", face = "bold")), 
    anc_plot + theme(axis.text.y = element_text(), axis.ticks.y = element_line()), 
    nrow = 1)

png("figures/GESMD_IWS_clustering/subgroup_exploration/clin_vars_subgroups_barplot_sup.png", width = 2000, height = 2000, res = 300)
clin_panel_sup
dev.off()



# clin_sum_plot <- clin_joint %>% 
#     select(all_of(clin_vars), sub_group, dataset) %>%
#     pivot_longer(cols = all_of(clin_vars), names_to = "Clin_var", values_to = "Value") %>%
#     filter(Clin_var != "PB_BLAST") %>%
#     group_by(Clin_var) %>%
#     mutate(sd_total = sd(Value, na.rm = TRUE)) %>%
#     ungroup() %>%
#     group_by(Clin_var, sub_group) %>%
#     summarize(m = mean(Value, na.rm = TRUE),
#         sd_total = first(sd_total)) %>% 
#     ungroup() %>%
#     group_by(Clin_var) %>%
#     mutate(d = (m - m[ sub_group == "Low blasts"])/sd_total,
#     Clin_var = factor(Clin_var, levels = c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT"))) %>%
#     filter(sub_group != "Low blasts") %>%
#     ggplot(aes(x = sub_group, y = Clin_var, fill = d )) +
#     geom_tile() +
#     scale_fill_gradient2(
#     low = "#0033ff",  
#     high = "#ff0000",
#     midpoint = 0
#     ) +
#   ylab("Clinical") +
#   xlab("Sub-groups") +
#   scale_x_discrete(position = "top") +
#   scale_y_discrete(limits = rev) +
#   theme_bw() +
#   theme(legend.position = "none")
# png("figures/GESMD_IWS_clustering/subgroup_exploration/clin_vars_subgroups_summary.png", width = 1800, height = 900, res = 300)
# clin_sum_plot
# dev.off()

## AGE

age_plot <-  ggplot(clin_joint3, aes(x = dataset, y = AGE, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("Age (Years)") +
      xlab("") +
      ggtitle("Age") +
      coord_cartesian(ylim = c(35, 100)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))


summary(lm(AGE ~ sub_group + dataset, mutate(clin_joint, sub_group = relevel(sub_group, "MDS-LB"))))
summary(lm(AGE ~ sub_group + dataset, mutate(clin_joint, sub_group = relevel(sub_group, "-7"))))
summary(lm(AGE ~ sub_group + dataset, mutate(clin_joint, sub_group = relevel(sub_group, "MDS-IB1"))))

## Sex
sex_plot <- clin_joint3 %>%
      group_by(sub_group, dataset) %>%
      summarize(FreqF = mean(SEX == "F")*100) %>%
      mutate(dataset = factor(dataset, levels = c("IWS", "GESMD", "MLL"))) %>%
      ggplot(aes(x = dataset, y = FreqF, fill = sub_group)) +
      geom_bar(stat = "identity") +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("Female Prop. (%)") +
      xlab("") +
      ggtitle("Sex") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
# png("figures/GESMD_IWS_clustering/subgroup_exploration/Sex_subgroups.png", width = 2000, height = 1000, res = 300)
# sex_plot
# dev.off()

png("figures/GESMD_IWS_clustering/subgroup_exploration/demographics_panel.png", width = 1500, height = 1600, res = 300)
plot_grid(age_plot + theme(legend.position = "none"), 
    sex_plot + theme(legend.position = "bottom"),
    ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1.3))
dev.off()

chisq.test(table(clin_joint$SEX, clin_joint$sub_group))
table(clin_joint$SEX, clin_joint$sub_group, clin_joint$dataset)
prop.table(table(clin_joint$SEX, clin_joint$sub_group, clin_joint$dataset), margin = c(2,3))




# sex_ORs <- lapply(sub_groups[sub_groups != "Low blasts"], function(group){
#         sel_df <- filter(clin_joint, sub_group %in% c(group, "Low blasts")) %>%
#             mutate(sub_group = factor(sub_group, levels = c("Low blasts", group)),
#             SEX = factor(SEX, levels = c("M", "F")))
#         tab <- table(sel_df$SEX, sel_df$sub_group)
#         fisher_res <- fisher.test(tab)
#         tibble(sub_group = group,
#                OR = fisher_res$estimate,
#                P_value = fisher_res$p.value)
#     }) %>%    Reduce(., f = rbind)



# dem_sum_plot <- bind_rows(clin_joint %>%
#     select(AGE, sub_group) %>%
#     mutate(sd_age = sd(AGE, na.rm = TRUE)) %>%
#     group_by(sub_group) %>%
#     summarize(m = mean(AGE, na.rm = TRUE),
#         sd = first(sd_age)) %>%
#     mutate(d = (m - m[sub_group == "Low blasts"])/sd,
#     var = "Age (Years)"),
# sex_ORs %>%
#     mutate(d = log(OR),
#     var = "Sex")
# ) %>% 
#     filter(sub_group != "Low blasts") %>%
#     mutate(sub_group = factor(sub_group, levels = sub_groups[sub_groups != "Low blasts"])) %>%
#     ggplot(aes(x = sub_group, y = var, fill = d )) +
#     geom_tile() +
#     scale_fill_gradient2(
#     low = "#0033ff",  
#     high = "#ff0000",
#     midpoint = 0
#     ) +
#   ylab("Demographics") +
#   xlab("Sub-groups") +
#   scale_x_discrete(position = "top") +
#   scale_y_discrete(limits = rev) +
#   theme_bw() +
#   theme(legend.position = "none")

# png("figures/GESMD_IWS_clustering/subgroup_exploration/demographics_subgroups_summary.png", width = 1800, height = 600, res = 300)
# dem_sum_plot
# dev.off()


## Mutations
# Select mutations Freq > 4% and not included in sub-groups
mutations <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "PHF6", "plus8", "delY", "del20q")

# plot_gene_frequency <- function(gene, scale_fill, output_dir = "figures/GESMD_IWS_clustering/subgroup_exploration/mutations/") {
#     png_filename <- file.path(output_dir, paste0(gene, "_subgroups.png"))
#    rbind(gesmd_dataset %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "GESMD"),
#       IWS_mds %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "IWS")) %>%
#       group_by(sub_group, dataset) %>%
#       summarize(Freq = mean(!!sym(gene))*100) %>%
#       mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
#      ggplot(aes(x = dataset, y = Freq, fill = sub_group)) +
#         geom_bar(stat = "identity") +
#         facet_grid(~ sub_group) +
#         scale_fill +
#         scale_y_continuous(name = "Frequency (%)", limits = c(0, 100)) +
#         xlab("") +
#         ggtitle(paste(gene, "frequency")) +
#         theme_bw() +
#         theme(
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#             strip.text = element_blank(), strip.background = element_blank(),
#             plot.title = element_text(hjust = 0.5)
#         )
#     ggsave(png_filename, width = 2000, height = 1000, dpi = 300, units = "px")
# }
# lapply(mutations, plot_gene_frequency, scale_fill = scale_fill)

# plot_gene_frequency2 <- function(gene, scale_fill) {
#     joint_df <- rbind(gesmd_dataset %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "GESMD"),
#       IWS_mds %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "IWS")) 

#     if (gene == "complex"){
#         joint_df[[gene]] <- ifelse(joint_df[[gene]] %in% c("complex", "1"), 1, 0)
#     }

#   df_plot <- joint_df %>%
#       group_by(sub_group, dataset) %>%
#       summarize(Freq = mean(!!sym(gene))*100) %>%
#       mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) 
#      if (gene == "TET2other") {
#          gene <- "TET2 Mono-allelic"
#      }
#      ggplot(df_plot, aes(x = dataset, y = Freq, fill = sub_group)) +
#         geom_bar(stat = "identity") +
#         facet_grid(~ sub_group) +
#         scale_fill +
#         scale_y_continuous(name = "Freq. (%)", limits = c(0, 100)) +
#         xlab("") +
#         ggtitle(gene) +
#         theme_bw() +
#         theme(
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#             strip.text = element_blank(), strip.background = element_blank(),
#             plot.title = element_text(hjust = 0.5),
#             legend.position = "none"
#         )
# }
# mut_plot <- lapply(mutations, plot_gene_frequency2, scale_fill = scale_fill)
# png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_panel.png", width = 4000, height = 2300, res = 300)
# plot_grid(
#     plot_grid(plotlist = mut_plot, ncol = 4, labels = "AUTO"),
#     legend, 
#     ncol = 1, rel_heights = c(1, 0.1)
# )
# dev.off()




joint_mut3 <- bind_rows(
    IWS_mds %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "GESMD"),
    hersh_mds %>% mutate(TET2other = TET2mono) %>% select(any_of(mutations), sub_group) %>% mutate(dataset = "MLL")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD", "MLL")))

joint_mut <- filter(joint_mut3, dataset != "MLL") %>%
    mutate(sub_group = droplevels(sub_group),
    dataset = droplevels(dataset))

hersh_mds$TET2other <- hersh_mds$TET2mono
mut_ORs <- lapply(mutations, function(var){
    lapply(seq_len(length(sub_groups) ), function(i){
        ref_cl <- sub_groups[i]
        lapply(sub_groups[sub_groups != ref_cl], function(comp_cl){
            sel_df <- filter(joint_mut, sub_group %in% c(ref_cl, comp_cl)) %>%
                mutate(sub_group = factor(sub_group, levels = c(ref_cl, comp_cl)))
            tab <- table(sel_df[[var]], sel_df$sub_group)
            fisher_res <- fisher.test(tab)

            sel_df_mll <- filter(hersh_mds, sub_group %in% c(ref_cl, comp_cl)) %>%
                mutate(sub_group = factor(sub_group, levels = c(ref_cl, comp_cl)))
            tab_mll <- table(sel_df_mll[[var]], sel_df_mll$sub_group)
            fisher_res_mll <- tryCatch(fisher.test(tab_mll), error = function(e) NULL)

            tibble(Gene = var,
                   Ref_cluster = ref_cl,
                   Comp_clust = comp_cl,
                   OR = fisher_res$estimate,
                   P_value = fisher_res$p.value,
                   OR_MLL = ifelse(is.null(fisher_res_mll), NA, fisher_res_mll$estimate),
                   P_value_MLL = ifelse(is.null(fisher_res_mll), NA, fisher_res_mll$p.value))
        }) %>% Reduce(., f = rbind)
    }) %>% Reduce(., f = rbind)
}) %>% Reduce(., f = rbind) %>%
    filter(OR > 1) %>%
    arrange(Gene, P_value)



write.table(mut_ORs, 
            file = "results/GESMD_IWS_clustering/subgroup_mutation_tests.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

mut_gene_plot <- c(
  "ASXL1", "SRSF2", "DNMT3A", "TET2other", 
  "RUNX1", "U2AF1", "BCOR", "ZRSR2", 
  "IDH2", "SETBP1", "DDX41", "PHF6",
  "plus8", "delY", "del20q"
)

mut_tib <- joint_mut3 %>%
    filter(sub_group %in% sub_groups) %>%
    mutate(sub_group = droplevels(sub_group)) %>%
    group_by(sub_group, dataset) %>%
    summarize_at(vars(all_of(mut_gene_plot)), mean, na.rm = TRUE) 

mut_mat <- as.matrix(mut_tib[, -c(1,2)])
rownames(mut_mat) <- seq_len(nrow(mut_tib))
colnames(mut_mat)[colnames(mut_mat) == "TET2other"] <- "TET2-mono"
colnames(mut_mat)[colnames(mut_mat) == "plus8"] <- "+8"
colnames(mut_mat)[colnames(mut_mat) == "delY"] <- "-Y"

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
          col = c(rep("black", length(unique(row_annot$sub_group)) - 1), "white") ,
          fill = colors
        ),        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", mut_mat[i, j] * 100), x, y, 
                    gp = gpar(fontsize = 10))
        },
        border = TRUE
)
png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_pheatmap.png", width = 2000, height = 2000, res = 300)
mut_heatmap
dev.off()

# mut_sum_plot <- mut_ORs %>%
#     filter(Comp_clust != "Low blasts") %>%
# #    filter(Gene %in% mutations[!mutations %in% c("IDH1", "TET2other")]) %>%
#     mutate(Gene = ifelse(Gene == "TET2other", "TET2 Mono-allelic", Gene),
#         Mut = factor(Gene, levels = c("RUNX1", "ASXL1", "ZRSR2", "SETBP1", "CBL", 
#     "IDH2", "BCOR", "TET2 Mono-allelic", "U2AF1", "IDH1", "DNMT3A", "SRSF2", "DDX41")), 
#     sub_group = factor(Comp_clust, levels = sub_groups)) %>%
#     mutate(OR = ifelse(OR == 0, 0.09, OR)) %>%
#     ggplot(aes(x = sub_group, y = Mut, fill = log(OR) )) +
#     geom_tile() +
#     scale_fill_gradient2(
#     low = "#0033ff",   # Azul para valores bajos
#    # mid = "white",     # Blanco para valores cercanos a la referencia
#     high = "#ff0000",  # Rojo para valores altos
#     midpoint = 0
#     ) +
#   ylab("Mutations") +
#   xlab("Sub-groups") +
#   scale_x_discrete(position = "top") +
#   scale_y_discrete(limits = rev) +
#   theme_bw() +
#   theme(legend.position = "none")

# png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_summary.png", width = 2000, height = 1200, res = 300)
# mut_sum_plot
# dev.off()

# kar_events <- list("plus8" = "8+", "delY" = "Y-", "del20q" = "del20q")

# joint_kar <- rbind(
#     IWS_mds %>% select(all_of(names(kar_events)), sub_group, CYTO_IPSSR) %>% mutate(dataset = "IWS"), 
#     gesmd_dataset %>% select(all_of(names(kar_events)), sub_group, CYTO_IPSSR) %>% mutate(dataset = "GESMD")
# ) %>%
#     mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")),
#     CYTO_IPSSR = fct_collapse(CYTO_IPSSR,
#     "Very-Good" = c("Very-Good", "Very Good"),
#     "Intermediate" = c("Intermediate", "Int"),
#     "Very-Poor" = c("Very-Poor", "Very Poor")))
# kar_ORs <- lapply(names(kar_events), function(var){
#     lapply(sub_groups[sub_groups != "MDS-LB"], function(group){
#         sel_df <- filter(joint_kar, sub_group %in% c(group, "MDS-LB")) %>%
#             mutate(sub_group = factor(sub_group, levels = c("MDS-LB", group)))
#         tab <- table(sel_df[[var]], sel_df$sub_group)
#         fisher_res <- fisher.test(tab)
#         tibble(Event = kar_events[[var]],
#                Sub_group = group,
#                OR = fisher_res$estimate,
#                P_value = fisher_res$p.value)
#     }) %>%    Reduce(., f = rbind)
# }) %>% Reduce(., f = rbind)

# write.table(kar_ORs, 
#             file = "results/GESMD_IWS_clustering/subgroup_karotype_tests.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             row.names = FALSE)



# kar_tib <- joint_mut3 %>%
#     filter(sub_group %in% sub_groups) %>%
#     mutate(sub_group = droplevels(sub_group)) %>%
#     group_by(sub_group, dataset) %>%
#     summarize_at(vars(all_of(names(kar_events))), mean, na.rm = TRUE) 

# mut_mat <- as.matrix(mut_tib[, -c(1,2)])
# rownames(mut_mat) <- seq_len(nrow(mut_tib))
# colnames(mut_mat)[colnames(mut_mat) == "TET2other"] <- "TET2-mono"

# row_annot <- data.frame(mut_tib[, 1])
# png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_pheatmap.png", width = 2000, height = 2000, res = 300)
# Heatmap(mut_mat, 
#         row_labels = mut_tib$dataset,
#         show_heatmap_legend = FALSE,
#         show_row_names = TRUE, 
#         show_column_names = TRUE,     
#         row_names_side = "left",
#         column_names_side = "top",
#         cluster_rows = FALSE, 
#         cluster_columns = TRUE, 
#         show_column_dend = FALSE, 
#         row_split = row_annot$sub_group, 
#         row_title = paste0("\n", unique(row_annot$sub_group)), 
#         row_gap = unit(2, "mm"), 
#         row_title_gp = gpar(
#           col = c(rep("black", length(unique(row_annot$sub_group)) - 1), "white") ,
#           fill = colors
#         ),        
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.1f", mut_mat[i, j] * 100), x, y, 
#                     gp = gpar(fontsize = 10))
#         },
#         border = TRUE
# )
# dev.off()



# kar_sum_plot <- kar_ORs %>%
#     filter(Sub_group != "Low blasts") %>%
#     mutate(Event = factor(Event, levels = c("8+", "Y-", "del20q")), 
#     sub_group = factor(Sub_group, levels = sub_groups),
#     OR = ifelse(OR == 0, 0.15, OR)) %>%
#     ggplot(aes(x = sub_group, y = Event, fill = log(OR) )) +
#     geom_tile() +
#     scale_fill_gradient2(
#     low = "#0033ff",   
#     high = "#ff0000", 
#     midpoint = 0
#     ) +
#   ylab("Karyotype") +
#   xlab("Sub-groups") +
#   scale_x_discrete(position = "top") +
#   scale_y_discrete(limits = rev) +
#   theme_bw() +
#   theme(legend.position = "none")

# png("figures/GESMD_IWS_clustering/subgroup_exploration/karyotype_subgroups_summary.png", width = 1800, height = 600, res = 300)
# kar_sum_plot
# dev.off()

# kar_plot <- lapply(names(kar_events), plot_gene_frequency2, scale_fill = scale_fill)
# kar_plot[[1]] <- kar_plot[[1]] + labs(title = "8+")
# kar_plot[[2]] <- kar_plot[[2]] + labs(title = "Y-")
# png("figures/GESMD_IWS_clustering/subgroup_exploration/karyotype_subgroups_panel.png", width = 2200, height = 1500, res = 300)
# plot_grid(
#     plot_grid(plotlist = kar_plot, ncol = 2, labels = "AUTO"),
#     legend, 
#     ncol = 1, rel_heights = c(1, 0.3)
# )
# dev.off()



joint_ipssm <- rbind(
    IWS_mds %>% select(IPSSM, CYTO_IPSSR, sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(IPSSM, CYTO_IPSSR, sub_group) %>% mutate(dataset = "GESMD"),
    hersh_mds %>% select(IPSSM, CYTO_IPSSR, sub_group) %>% 
    mutate(dataset = "MLL", IPSSM = gsub(" ", "-", IPSSM))
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD", "MLL"))) 
ipssm_summary <- joint_ipssm %>% filter(!is.na(IPSSM) & sub_group %in% sub_groups) %>%
  mutate(sub_group = droplevels(sub_group)) %>%
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
  ggplot(aes(y = fct_rev(dataset), x = Freq*100, fill = fct_rev(IPSSM))) +
  geom_bar(stat = "identity", orientation = "y") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "IPSSM classification",
    y = "",
    x = "",
    fill = "IPSSM") +
         facet_grid2(sub_group ~ ., 
                    switch = "y", 
                    scales = "free_y",
                    independent = "y",
                    strip = strip_themed(
                        background_y = elem_list_rect(fill = colors)
                )) +
  theme(plot.title = element_text(hjust = 0.5),
  legend.position = "top", strip.background = element_rect(),
  strip.text.y.left = element_text(angle = 0, color = "white", face = "bold")) +
  scale_y_discrete(position = "right")
png("figures/GESMD_IWS_clustering/subgroup_exploration/IPSSM_summary.png", res = 300, height = 900, width = 3000)
ipssm_summary
dev.off()


ht_grob <- grid.grabExpr(draw(mut_heatmap))
png("figures/GESMD_IWS_clustering/subgroup_exploration/panel_summary.png", res = 300, height = 4300, width = 3500)
plot_grid(
    plot_grid(clin_panel, ipssm_summary, ncol = 2, rel_widths = c(2, 1), labels = c("A", "C")),
    ht_grob, ncol = 1, labels = c("", "B"), rel_heights = c(1, 1)
)
dev.off()

ipssrcyto_summary <- joint_ipssm %>%
    filter(!is.na(CYTO_IPSSR)  & sub_group %in% sub_groups) %>%
    mutate(CYTO_IPSSR = fct_recode(CYTO_IPSSR,
    "Very Good" = "Very-Good",
    "Int" = "Intermediate",
    "Very Poor" = "Very-Poor")) %>%
  group_by(sub_group, CYTO_IPSSR, dataset) %>%
  summarize(N = n()) %>%
  group_by(sub_group, dataset) %>%
  mutate(Freq = N/sum(N),
  dataset = factor(dataset, levels = c("MLL", "GESMD", "IWS"))) %>%
  ggplot(aes(x = dataset, y = Freq*100, fill = CYTO_IPSSR)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fee08b", "#66bd63", "#2ca25f")) +
  labs(
    title = "Cytogenetic Risk category",
    y = "Individuals (%)",
    x = "",
    fill = "") +
    facet_grid(sub_group ~ .) +
  theme(plot.title = element_text(hjust = 0.5),
  legend.position = "bottom",  strip.background = element_rect()) +
  coord_flip() 
png("figures/GESMD_IWS_clustering/subgroup_exploration/IPSSR_CYTO_summary.png", res = 300, height = 1800, width = 1500)
ipssrcyto_summary
dev.off()


joint_ipssm2 <- subset(joint_ipssm, dataset != "MLL" & sub_group %in% sub_groups) %>%
    mutate(sub_group = droplevels(sub_group),
    dataset = droplevels(dataset)) %>%
    mutate(CYTO_IPSSR = fct_recode(CYTO_IPSSR,
    "Very Good" = "Very-Good",
    "Intermediate" = "Int",
    "Very Poor" = "Very-Poor"))
table(joint_ipssm2$CYTO_IPSSR, joint_ipssm2$sub_group)
chisq.test(table(joint_ipssm2$CYTO_IPSSR, joint_ipssm2$sub_group))

joint_ipssm2_f <- filter(joint_ipssm2, sub_group != "-7") %>%
    mutate(sub_group = droplevels(sub_group))
table(joint_ipssm2_f$CYTO_IPSSR %in% c("Very Good", "Good"), joint_ipssm2_f$sub_group)
chisq.test(table(joint_ipssm2_f$CYTO_IPSSR %in% c("Very Good", "Good"), joint_ipssm2_f$sub_group))


lapply(sub_groups[sub_groups != "-7"], function(group){
    print(group)
    print(fisher.test(table(!joint_ipssm2_f$CYTO_IPSSR %in% c("Very Good", "Good"), 
        joint_ipssm2_f$sub_group == group)))
})

joint_ipssm2_m <- filter(joint_ipssm2, sub_group %in% c("MDS-LB", "MDS-IB1", "MDS-IB2")) %>%
    mutate(sub_group = droplevels(sub_group))

fisher.test(table(!joint_ipssm2_f$CYTO_IPSSR %in% c("Very Good", "Good"), 
        joint_ipssm2_f$sub_group == "MDS-LB"))


# png("figures/GESMD_IWS_clustering/subgroup_exploration/overall_summary.png", res = 300, height = 2200, width = 3800)
# plot_grid(
#     plot_grid(
#         plot_grid(clin_sum_plot, kar_sum_plot,
#             ncol = 1, labels = c("A", "B"), rel_heights = c(1.5, 1)),
#         mut_sum_plot, ncol = 2, labels = c("", "D")),
#     ipssrcyto_summary,
#     ipssm_summary, 
#     ncol = 1, labels = c("",  "C", "E"), rel_heights = c(3, 2, 2)
# )
# dev.off()


# ## Make descriptives for each dataset per cluster
# getIQR <- function(vec){
#    sprintf("%.1f (%.1f-%.1f)", 
#            median(vec, na.rm = TRUE), 
#            quantile(vec, probs = 0.25, na.rm = TRUE),
#            quantile(vec, probs = 0.75, na.rm = TRUE))
# }


# getProp <- function(cond, var){
#     sprintf("%i (%.1f%%)", sum(cond, na.rm = TRUE), mean(!is.na(var) & cond)*100)
# }

# getPropNA <- function(var){
#     sprintf("%i (%.1f%%)", sum(is.na(var)), mean(is.na(var))*100)
# }
# summarize_fun <- function(df){
#      summarize(df, 
#             N_dataset = first(N),
#             N = sprintf("%i (%.1f%%)", n(), n()/N_dataset*100),
#             Females = getProp(SEX == "F", SEX),
#             Males = getProp(SEX == "M", SEX),
#             Age = getIQR(AGE),
#             `BM Blasts` = getIQR(BM_BLAST),
#             `WBC count` = getIQR(WBC),
#             `Neutrophil Count` = getIQR(ANC),
#             `Monocyte Count` = getIQR(MONOCYTES),
#             HB = getIQR(HB),
#             PLT = getIQR(PLT),
#             `Low blasts` = getProp(consensus == "Low blasts", consensus),
#             `MDS-IB1` = getProp(consensus == "MDS-IB1", consensus),
#             `MDS-IB2` = getProp(consensus == "MDS-IB2", consensus),
#             `Very-Low` = getProp(IPSSM == "Very-Low", IPSSM),
#             `Low` = getProp(IPSSM == "Low", IPSSM),
#             `Moderate-Low` = getProp(IPSSM == "Moderate-Low", IPSSM),
#             `Moderate-High` = getProp(IPSSM == "Moderate-High", IPSSM),
#             `High` = getProp(IPSSM == "High", IPSSM),
#             `Very-High` = getProp(IPSSM == "Very-High",  IPSSM),
#             IPSSM_NA = getPropNA(IPSSM),
#             `8+` = getProp(plus8 == 1, plus8),
#             `Y-` = getProp(delY == 1, delY),
#             del20q = getProp(del20q == 1, del20q),
#             ASXL1 = getProp(ASXL1 == 1, ASXL1),
#             SRSF2 = getProp(SRSF2 == 1, SRSF2),
#             DNMT3A = getProp(DNMT3A == 1, DNMT3A),
#             RUNX1 = getProp(RUNX1 == 1, RUNX1),
#             U2AF1 = getProp(U2AF1 == 1, U2AF1),
#             BCOR = getProp(BCOR == 1, BCOR),
#             ZRSR2 = getProp(ZRSR2 == 1, ZRSR2),
#             TET2other = getProp(TET2other == 1, TET2other),
#             IDH1 = getProp(IDH1 == 1, IDH1),
#             IDH2 = getProp(IDH2 == 1, IDH2),
#             SETBP1 = getProp(SETBP1 == 1, SETBP1),
#             CBL = getProp(CBL == 1, CBL),
#             DDX41 = getProp(DDX41 == 1, DDX41)    
#            )
# }

# ### Combined
# joint_subgroup_descriptives <- bind_rows(gesmd_dataset %>% mutate(dataset = "GESMD"),
#           IWS_mds %>% mutate(dataset = "IWS", complex = ifelse(complex == "complex", 1, 0))) %>%
#     mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
#     group_by(dataset) %>%
#     mutate(N = n()) %>%
#     group_by(dataset, sub_group) %>%
#     summarize_fun() %>%
#     select(-N_dataset) %>%
#     arrange(sub_group, dataset) %>%
#     t() 

# write.table(joint_subgroup_descriptives, 
#             file = "results/GESMD_IWS_clustering/joint_subgroup_descriptives.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             col.names = NA)


# ## Tests
# poisson_test <- function(var, df){
#     poiss_base <- glm(formula (paste(var, " ~ 1")), df, 
#                       family = "poisson")
#     poiss_lm <- glm(formula (paste(var, " ~ sub_group")), df, 
#                       family = "poisson")
#     anova_res <- anova(poiss_base, poiss_lm, test = "Chisq")                      
#     anova_res$`Pr(>Chi)`[2]
# }

# lm_test <- function(var, df){
#     lm_res <- summary(lm(formula (paste(var, " ~ sub_group")), df))
#     fstats <- lm_res$fstatistic
#     pf(fstats[1], fstats[2], fstats[3], lower.tail = FALSE)
#  }

# chisq_test <- function(var, df){
#     chisq_res <- chisq.test(table(df[[var]], df$sub_group))
#     chisq_res$p.value
# }

# # Perform tests
# gesmd_test <- c(sapply(c("SEX", "consensus", "IPSSM", mutations), chisq_test, df = gesmd_dataset),
#             sapply(c("AGE", "HB", "PLT"), lm_test, df = gesmd_dataset),
#             sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = gesmd_dataset))
# names(gesmd_test) <- c("Females", "Low blasts", "IPSSM_NA", mutations, "Age", "HB", "PLT", 
#     "BM Blasts", "WBC count", "Neutrophil Count", "Monocyte Count") 

# gesmd_subgroup_df$p_value <- gesmd_test[rownames(gesmd_subgroup_df)]

# write.table(gesmd_subgroup_df, 
#             file = "results/GESMD_IWS_clustering/gesmd_subgroup_descriptives.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             col.names = NA)

# ### IWS
# iws_subgroup <- rbind(IWS_mds  %>%
#     group_by(sub_group) %>%
#     summarize_fun(),
#     cbind(tibble(sub_group = "Total"),  summarize_fun(IWS_mds))
# ) %>%
#     t() 
# iws_subgroup_df <- as.data.frame(iws_subgroup[-1, ]) 
# colnames(iws_subgroup_df) <- iws_subgroup[1, ]
# rownames(iws_subgroup_df) <- rownames(iws_subgroup)[-1]



# # Perform tests
# iws_test <- c(sapply(c("SEX", "consensus", "IPSSM", mutations), chisq_test, df = IWS_mds),
#             sapply(c("AGE", "HB", "PLT"), lm_test, df = IWS_mds),
#             sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = IWS_mds))
# names(iws_test) <- c("Females", "Low blasts", "IPSSM_NA", mutations, "Age", "HB", "PLT", 
#     "BM Blasts", "WBC count", "Neutrophil Count", "Monocyte Count") 

# iws_subgroup_df$p_value <- iws_test[rownames(iws_subgroup_df)]

# write.table(iws_subgroup_df, 
#             file = "results/GESMD_IWS_clustering/IWS_subgroup_descriptives.txt", 
#             sep = "\t", 
#             quote = FALSE, 
#             col.names = NA)
