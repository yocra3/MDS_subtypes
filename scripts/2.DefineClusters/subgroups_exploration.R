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
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(MASS)
library(tidyverse)
library(cowplot)

load("results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

# colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
#     "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")
colors6 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#999999", "grey40",  "black")

scale_col <- scale_color_manual(values = colors6, name = "Sub-group")
scale_fill <- scale_fill_manual(values = colors6, name = "Sub-group")

## Proportion of each subgroup in each dataset
mean(IWS_mds$sub_group %in% c("EZH2", "TET2 bi-allelic", "7-", "STAG2"))
mean(gesmd_dataset$sub_group %in% c("EZH2", "TET2 bi-allelic", "7-", "STAG2"))

prop.table(table(IWS_mds$sub_group))
prop.table(table(gesmd_dataset$sub_group))

## Clinical variables
clusters <- levels(IWS_mds$sub_group)

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


### WBC
png("figures/GESMD_IWS_clustering/subgroup_exploration/WBC_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("WBC", "WBC count")
dev.off()


### BM_BLAST
png("figures/GESMD_IWS_clustering/subgroup_exploration/BM_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("BM_BLAST", "BM Blasts proportion")
dev.off()

### ANC
png("figures/GESMD_IWS_clustering/subgroup_exploration/ANC_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("ANC", "ANC")
dev.off()

### MONOCYTES
png("figures/GESMD_IWS_clustering/subgroup_exploration/MONOCYTES_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("MONOCYTES", "Monocytes")
dev.off()

### PLT
png("figures/GESMD_IWS_clustering/subgroup_exploration/PLT_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("PLT", "Platelets")
dev.off()

### PB_BLAST
png("figures/GESMD_IWS_clustering/subgroup_exploration/PB_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("PB_BLAST", "PB Blasts proportion")
dev.off()

### HB
png("figures/GESMD_IWS_clustering/subgroup_exploration/HB_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("HB", "Hemoglobin")
dev.off()


png("figures/GESMD_IWS_clustering/subgroup_exploration/WBC_subgroups.png", width = 2000, height = 1200, res = 300)
makeBoxPlots("WBC", "WBC count")
dev.off()


### Panel of clinical variables
clin_names <- list(BM_BLAST = "BM Blasts",
                   WBC = "WBC count",
                   ANC = "ANC",
                   MONOCYTES = "Monocytes",
                   PLT = "Platelets",
                   HB = "Hemoglobin")
clin_plots <- lapply(names(clin_names), function(var){
    makeBoxPlots(var, clin_names[[var]]) +
    theme(legend.position = "none")
})
clin_plots[[1]] <- clin_plots[[1]] + labs(y = "BM blasts (%)")

legend <- get_plot_component(makeBoxPlots("ANC", "ANC") + 
theme(legend.position = "bottom"),
 "guide-box", return_all = TRUE)[[3]]

png("figures/GESMD_IWS_clustering/subgroup_exploration/clin_vars_panel.png", width = 2000, height = 1900, res = 300)
plot_grid(
    plot_grid(plotlist = clin_plots, ncol = 2,
        labels = "AUTO"),
    legend,
    ncol = 1, rel_heights = c(1, 0.1))
dev.off()


## Tests
sub_groups <- c("EZH2", "TET2 bi-allelic", "7-", "STAG2", "Low blasts", "MDS-IB1", "MDS-IB2")
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_joint <- bind_rows(
    IWS_mds %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(clin_vars), sub_group, AGE, SEX) %>% mutate(dataset = "GESMD")
)

nb <- c("BM_BLAST", "PB_BLAST")
res_nb <- lapply(nb, function(var){
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm.nb(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp))
     coefs <- poiss_lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
              filter(Comp_clust %in% sub_groups)
     colnames(coef_df)[2] <- "P_value"
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value)
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
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp, 
                      family = "poisson"))
     coefs <- poiss_lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate)) %>%
              filter(Comp_clust %in% sub_groups)
     colnames(coef_df)[2] <- "P_value"
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value)
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
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clin_joint, cluster = relevel(sub_group, cl))
     lm <- summary(lm(formula (paste(var, " ~ cluster + AGE + SEX + dataset")), tmp))
     coefs <- lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = Estimate) %>%
              filter(Comp_clust %in% sub_groups)
     colnames(coef_df)[2] <- "P_value"
     
     select(coef_df, Ref_cluster, Comp_clust, Stat, P_value)
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
    mutate(d = (m - m[ sub_group == "Low blasts"])/sd_total,
    Clin_var = factor(Clin_var, levels = c("BM_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT"))) %>%
    filter(sub_group != "Low blasts") %>%
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
  theme(legend.position = "none")
png("figures/GESMD_IWS_clustering/subgroup_exploration/clin_vars_subgroups_summary.png", width = 1800, height = 900, res = 300)
clin_sum_plot
dev.off()

## AGE
age_plot <- makeBoxPlots("AGE", "Age") + 
labs(y = "Age (Years)")


png("figures/GESMD_IWS_clustering/subgroup_exploration/AGE_subgroups.png", width = 2000, height = 1200, res = 300)
age_plot
dev.off()


summary(lm(AGE ~ sub_group + dataset, mutate(clin_joint, sub_group = relevel(sub_group, "Low blasts"))))

## Sex
sex_plot <- clin_joint %>%
      group_by(sub_group, dataset) %>%
      summarize(FreqF = mean(SEX == "F")*100) %>%
      mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
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
png("figures/GESMD_IWS_clustering/subgroup_exploration/Sex_subgroups.png", width = 2000, height = 1000, res = 300)
sex_plot
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_exploration/demographics_panel.png", width = 1500, height = 1600, res = 300)
plot_grid(age_plot + theme(legend.position = "none"), 
    sex_plot + theme(legend.position = "none"),
    legend,
    ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1, 0.3))
dev.off()

chisq.test(table(clin_joint$SEX, clin_joint$sub_group))
table(clin_joint$SEX, clin_joint$sub_group, clin_joint$dataset)
prop.table(table(clin_joint$SEX, clin_joint$sub_group, clin_joint$dataset), margin = c(2,3))




sex_ORs <- lapply(sub_groups[sub_groups != "Low blasts"], function(group){
        sel_df <- filter(clin_joint, sub_group %in% c(group, "Low blasts")) %>%
            mutate(sub_group = factor(sub_group, levels = c("Low blasts", group)),
            SEX = factor(SEX, levels = c("M", "F")))
        tab <- table(sel_df$SEX, sel_df$sub_group)
        fisher_res <- fisher.test(tab)
        tibble(sub_group = group,
               OR = fisher_res$estimate,
               P_value = fisher_res$p.value)
    }) %>%    Reduce(., f = rbind)



dem_sum_plot <- bind_rows(clin_joint %>%
    select(AGE, sub_group) %>%
    mutate(sd_age = sd(AGE, na.rm = TRUE)) %>%
    group_by(sub_group) %>%
    summarize(m = mean(AGE, na.rm = TRUE),
        sd = first(sd_age)) %>%
    mutate(d = (m - m[sub_group == "Low blasts"])/sd,
    var = "Age (Years)"),
sex_ORs %>%
    mutate(d = log(OR),
    var = "Sex")
) %>% 
    filter(sub_group != "Low blasts") %>%
    mutate(sub_group = factor(sub_group, levels = sub_groups[sub_groups != "Low blasts"])) %>%
    ggplot(aes(x = sub_group, y = var, fill = d )) +
    geom_tile() +
    scale_fill_gradient2(
    low = "#0033ff",  
    high = "#ff0000",
    midpoint = 0
    ) +
  ylab("Demographics") +
  xlab("Sub-groups") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_exploration/demographics_subgroups_summary.png", width = 1800, height = 600, res = 300)
dem_sum_plot
dev.off()


## Mutations
# Select mutations Freq > 4% and not included in models
mutations <- c("ASXL1", "SRSF2", "DNMT3A", "TET2other", "RUNX1", "U2AF1",  
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "CBL", "IDH1")

plot_gene_frequency <- function(gene, scale_fill, output_dir = "figures/GESMD_IWS_clustering/subgroup_exploration/mutations/") {
    png_filename <- file.path(output_dir, paste0(gene, "_subgroups.png"))
   rbind(gesmd_dataset %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "GESMD"),
      IWS_mds %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "IWS")) %>%
      group_by(sub_group, dataset) %>%
      summarize(Freq = mean(!!sym(gene))*100) %>%
      mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
     ggplot(aes(x = dataset, y = Freq, fill = sub_group)) +
        geom_bar(stat = "identity") +
        facet_grid(~ sub_group) +
        scale_fill +
        scale_y_continuous(name = "Frequency (%)", limits = c(0, 100)) +
        xlab("") +
        ggtitle(paste(gene, "frequency")) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)
        )
    ggsave(png_filename, width = 2000, height = 1000, dpi = 300, units = "px")
}
lapply(mutations, plot_gene_frequency, scale_fill = scale_fill)

plot_gene_frequency2 <- function(gene, scale_fill) {
   rbind(gesmd_dataset %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "GESMD"),
      IWS_mds %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "IWS")) %>%
      group_by(sub_group, dataset) %>%
      summarize(Freq = mean(!!sym(gene))*100) %>%
      mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
     ggplot(aes(x = dataset, y = Freq, fill = sub_group)) +
        geom_bar(stat = "identity") +
        facet_grid(~ sub_group) +
        scale_fill +
        scale_y_continuous(name = "Freq. (%)", limits = c(0, 100)) +
        xlab("") +
        ggtitle(gene) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
        )
}
mut_plot <- lapply(mutations[!mutations %in% c("IDH1", "TET2other")], plot_gene_frequency2, scale_fill = scale_fill)
png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_panel.png", width = 2400, height = 2300, res = 300)
plot_grid(
    plot_grid(plotlist = mut_plot, ncol = 3, labels = "AUTO"),
    legend, 
    ncol = 1, rel_heights = c(1, 0.1)
)
dev.off()




joint_mut <- rbind(
    IWS_mds %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(mutations), sub_group) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")))

mut_ORs <- lapply(mutations, function(var){
    lapply(sub_groups[sub_groups != "Low blasts"], function(group){
        sel_df <- filter(joint_mut, sub_group %in% c(group, "Low blasts")) %>%
            mutate(sub_group = factor(sub_group, levels = c("Low blasts", group)))
        tab <- table(sel_df[[var]], sel_df$sub_group)
        fisher_res <- fisher.test(tab)
        tibble(Gene = var,
               Comp_clust = group,
               OR = fisher_res$estimate,
               P_value = fisher_res$p.value)
    }) %>%    Reduce(., f = rbind)
}) %>% Reduce(., f = rbind)

write.table(mut_ORs, 
            file = "results/GESMD_IWS_clustering/subgroup_mutation_tests.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)



mut_sum_plot <- mut_ORs %>%
    filter(Comp_clust != "Low blasts") %>%
    filter(Gene %in% mutations[!mutations %in% c("IDH1", "TET2other")]) %>%
    mutate(Mut = factor(Gene, levels = c("RUNX1", "ASXL1", "ZRSR2", "SETBP1", "CBL", 
    "IDH2", "BCOR", "U2AF1", "DNMT3A", "SRSF2", "DDX41")), 
    sub_group = factor(Comp_clust, levels = sub_groups)) %>%
    ggplot(aes(x = sub_group, y = Mut, fill = log(OR) )) +
    geom_tile() +
    scale_fill_gradient2(
    low = "#0033ff",   # Azul para valores bajos
   # mid = "white",     # Blanco para valores cercanos a la referencia
    high = "#ff0000",  # Rojo para valores altos
    midpoint = 0
    ) +
  ylab("Mutations") +
  xlab("Sub-groups") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_exploration/mutation_subgroups_summary.png", width = 1800, height = 1200, res = 300)
mut_sum_plot
dev.off()

kar_events <- list("plus8" = "8+", "delY" = "Y-", "del20q" = "del20q")

joint_kar <- rbind(
    IWS_mds %>% select(all_of(names(kar_events)), sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(all_of(names(kar_events)), sub_group) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD")))
kar_ORs <- lapply(names(kar_events), function(var){
    lapply(sub_groups[sub_groups != "Low blasts"], function(group){
        sel_df <- filter(joint_kar, sub_group %in% c(group, "Low blasts")) %>%
            mutate(sub_group = factor(sub_group, levels = c("Low blasts", group)))
        tab <- table(sel_df[[var]], sel_df$sub_group)
        fisher_res <- fisher.test(tab)
        tibble(Event = kar_events[[var]],
               Sub_group = group,
               OR = fisher_res$estimate,
               P_value = fisher_res$p.value)
    }) %>%    Reduce(., f = rbind)
}) %>% Reduce(., f = rbind)

write.table(kar_ORs, 
            file = "results/GESMD_IWS_clustering/subgroup_karotype_tests.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

kar_sum_plot <- kar_ORs %>%
    filter(Sub_group != "Low blasts") %>%
    mutate(Event = factor(Event, levels = c("8+", "Y-", "del20q")), 
    sub_group = factor(Sub_group, levels = sub_groups)) %>%
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
  theme(legend.position = "none")

png("figures/GESMD_IWS_clustering/subgroup_exploration/karyotype_subgroups_summary.png", width = 1800, height = 600, res = 300)
kar_sum_plot
dev.off()

kar_plot <- lapply(names(kar_events), plot_gene_frequency2, scale_fill = scale_fill)
kar_plot[[1]] <- kar_plot[[1]] + labs(title = "8+")
kar_plot[[2]] <- kar_plot[[2]] + labs(title = "Y-")
png("figures/GESMD_IWS_clustering/subgroup_exploration/karyotype_subgroups_panel.png", width = 2400, height = 700, res = 300)
plot_grid(
    plot_grid(plotlist = kar_plot, ncol = 3, labels = "AUTO"),
    legend, 
    ncol = 1, rel_heights = c(1, 0.3)
)
dev.off()



ipssm_summary <- rbind(
    IWS_mds %>% select(IPSSM, sub_group) %>% mutate(dataset = "IWS"), 
    gesmd_dataset %>% select(IPSSM, sub_group) %>% mutate(dataset = "GESMD")
) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
  filter(!is.na(IPSSM)) %>%
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
    x = "Prognostic group",
    fill = "IPSSM") +
    facet_grid(~ sub_group) +
  theme(plot.title = element_text(hjust = 0.5)) 
png("figures/GESMD_IWS_clustering/subgroup_exploration/IPSSM_summary.png", res = 300, height = 1200, width = 1800)
ipssm_summary
dev.off()

png("figures/GESMD_IWS_clustering/subgroup_exploration/overall_summary.png", res = 300, height = 2200, width = 2800)
plot_grid(
    plot_grid(
        plot_grid(clin_sum_plot, kar_sum_plot,
            ncol = 1, labels = c("A", "C"), rel_heights = c(2, 1)),
        plot_grid(dem_sum_plot, mut_sum_plot,
            ncol = 1, labels = c("B", "D"), rel_heights = c(1, 2)),
        ncol = 2),
    ipssm_summary, 
    ncol = 1, labels = c("", "E"), rel_heights = c(3, 2)
)
dev.off()


## Make descriptives for each dataset per cluster
getIQR <- function(vec){
   sprintf("%.1f (%.1f-%.1f)", 
           median(vec, na.rm = TRUE), 
           quantile(vec, probs = 0.25, na.rm = TRUE),
           quantile(vec, probs = 0.75, na.rm = TRUE))
}


getProp <- function(cond, var){
    sprintf("%i (%.1f%%)", sum(cond, na.rm = TRUE), mean(!is.na(var) & cond)*100)
}

getPropNA <- function(var){
    sprintf("%i (%.1f%%)", sum(is.na(var)), mean(is.na(var))*100)
}
summarize_fun <- function(df){
     summarize(df, 
            N_dataset = first(N),
            N = sprintf("%i (%.1f%%)", n(), n()/N_dataset*100),
            Females = getProp(SEX == "F", SEX),
            Males = getProp(SEX == "M", SEX),
            Age = getIQR(AGE),
            `BM Blasts` = getIQR(BM_BLAST),
            `WBC count` = getIQR(WBC),
            `Neutrophil Count` = getIQR(ANC),
            `Monocyte Count` = getIQR(MONOCYTES),
            HB = getIQR(HB),
            PLT = getIQR(PLT),
            `Low blasts` = getProp(consensus == "Low blasts", consensus),
            `MDS-IB1` = getProp(consensus == "MDS-IB1", consensus),
            `MDS-IB2` = getProp(consensus == "MDS-IB2", consensus),
            `Very-Low` = getProp(IPSSM == "Very-Low", IPSSM),
            `Low` = getProp(IPSSM == "Low", IPSSM),
            `Moderate-Low` = getProp(IPSSM == "Moderate-Low", IPSSM),
            `Moderate-High` = getProp(IPSSM == "Moderate-High", IPSSM),
            `High` = getProp(IPSSM == "High", IPSSM),
            `Very-High` = getProp(IPSSM == "Very-High",  IPSSM),
            IPSSM_NA = getPropNA(IPSSM),
            `8+` = getProp(plus8 == 1, plus8),
            `Y-` = getProp(delY == 1, delY),
            del20q = getProp(del20q == 1, del20q),
            ASXL1 = getProp(ASXL1 == 1, ASXL1),
            SRSF2 = getProp(SRSF2 == 1, SRSF2),
            DNMT3A = getProp(DNMT3A == 1, DNMT3A),
            RUNX1 = getProp(RUNX1 == 1, RUNX1),
            U2AF1 = getProp(U2AF1 == 1, U2AF1),
            BCOR = getProp(BCOR == 1, BCOR),
            ZRSR2 = getProp(ZRSR2 == 1, ZRSR2),
            IDH2 = getProp(IDH2 == 1, IDH2),
            SETBP1 = getProp(SETBP1 == 1, SETBP1),
            CBL = getProp(CBL == 1, CBL),
            DDX41 = getProp(DDX41 == 1, DDX41)    
           )
}

### Combined
joint_subgroup_descriptives <- bind_rows(gesmd_dataset %>% mutate(dataset = "GESMD") %>% select(-complex),
          IWS_mds %>% mutate(dataset = "IWS") %>% select(-complex)) %>%
    mutate(dataset = factor(dataset, levels = c("IWS", "GESMD"))) %>%
    group_by(dataset) %>%
    mutate(N = n()) %>%
    group_by(dataset, sub_group) %>%
    summarize_fun() %>%
    select(-N_dataset) %>%
    arrange(sub_group, dataset) %>%
    t() 

write.table(joint_subgroup_descriptives, 
            file = "results/GESMD_IWS_clustering/joint_subgroup_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)


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
