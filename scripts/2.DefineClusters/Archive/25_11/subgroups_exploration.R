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

load( "results/GESMD_IWS_clustering/gesmd_IWS_full.Rdata")

colors <- c("black", "grey40", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#0000FF")
scale_col <- scale_color_manual(values = colors, name = "Sub-group")
scale_fill <- scale_fill_manual(values = colors, name = "Sub-group")

change_levels <- function(df){
    df <- mutate(df, sub_group = factor(sub_group,
        levels = c("Highly Leukopenic", "Midly Leukopenic", "TET2 monoallelic", 
            "TET2 bi-allelic", "Y-", "8+", "7-", "del20q", "del7q", "complex", "STAG2")))
    return(df)
}

## Clinical variables
clinical_blasts <- change_levels(clinical_blasts)
clusters <- levels(clinical_blasts$sub_group)

### WBC
png("figures/GESMD_IWS_clustering/subgroup_exploration/WBC_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, WBC) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, WBC) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = WBC, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("WBC count") +
      xlab("") +
      ggtitle("WBC count") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()


### BM_BLAST
png("figures/GESMD_IWS_clustering/subgroup_exploration/BM_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, BM_BLAST) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, BM_BLAST) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = BM_BLAST, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("BM (%)") +
      xlab("") +
      ggtitle("BM Blasts proportion") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()

### ANC
png("figures/GESMD_IWS_clustering/subgroup_exploration/ANC_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, ANC) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, ANC) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = ANC, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("ANC") +
      xlab("") +
      ggtitle("ANC") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()
### MONOCYTES
png("figures/GESMD_IWS_clustering/subgroup_exploration/MONOCYTES_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, MONOCYTES) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, MONOCYTES) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = MONOCYTES, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("MONOCYTES") +
      xlab("") +
      ggtitle("MONOCYTES") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()
### PLT
png("figures/GESMD_IWS_clustering/subgroup_exploration/PLT_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, PLT) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, PLT) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = PLT, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("PLT") +
      xlab("") +
      ggtitle("PLT") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()
### PB_BLAST
png("figures/GESMD_IWS_clustering/subgroup_exploration/PB_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, PB_BLAST) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, PB_BLAST) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = PB_BLAST, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("PB (%)") +
      xlab("") +
      ggtitle("PB Blasts proportion") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()
### HB
png("figures/GESMD_IWS_clustering/subgroup_exploration/HB_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, HB) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, HB) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = HB, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("HB") +
      xlab("") +
      ggtitle("Hemoglobin") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()

## Tests
nb <- c("BM_BLAST", "PB_BLAST")
res_nb <- lapply(nb, function(var){
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clinical_blasts, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm.nb(formula (paste(var, " ~ cluster")), tmp))
     coefs <- poiss_lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate))
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

poisson <- c("WBC", "ANC", "MONOCYTES")
res_poisson <- lapply(poisson, function(var){
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clinical_blasts, cluster = relevel(sub_group, cl))
     poiss_lm <- summary(glm(formula (paste(var, " ~ cluster")), tmp, 
                      family = "poisson"))
     coefs <- poiss_lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = exp(Estimate))
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

normal <- c("HB", "PLT")
res_norm <- lapply(normal, function(var){
    cell_test <- lapply(seq_len(length(clusters) - 1), function(i){
      cl <- clusters[i]
     tmp <- mutate(clinical_blasts, cluster = relevel(sub_group, cl))
     lm <- summary(lm(formula (paste(var, " ~ cluster ")), tmp))
     coefs <- lm$coefficients[-1, ]
     coef_df <- as_tibble(coefs[, c(1, 4)]) %>%
       mutate(Ref_cluster = cl, 
              Comp_clust = gsub("cluster", "", rownames(coefs)),
              Stat = Estimate)
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

## AGE
png("figures/GESMD_IWS_clustering/subgroup_exploration/AGE_subgroups.png", width = 2000, height = 1200, res = 300)
rbind(gesmd %>% select(sub_group, AGE) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, AGE) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      ggplot(aes(x = dataset, y = AGE, fill = sub_group)) +
      geom_boxplot() +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("Age (Years)") +
      xlab("") +
      ggtitle("Age") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()
summary(lm(AGE ~ sub_group, clinical_blasts))

## Sex
png("figures/GESMD_IWS_clustering/subgroup_exploration/Sex_subgroups.png", width = 2000, height = 1000, res = 300)
rbind(gesmd %>% select(sub_group, SEX) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, SEX) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      group_by(sub_group, dataset) %>%
      summarize(FreqF = mean(SEX == "F")*100) %>%
      ggplot(aes(x = dataset, y = FreqF, fill = sub_group)) +
      geom_bar(stat = "identity") +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("Female proportion (%)") +
      xlab("") +
      ggtitle("Sex differences") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()

clin_noY <- subset(clinical_blasts, sub_group != "Y-") %>%
    mutate(sub_group = droplevels(sub_group))
chisq.test(table(clin_noY$SEX, clin_noY$sub_group))

## Mutations
# Select: ASXL1, SRSF2, DNMT3A, RUNX1, U2AF1, EZH2, BCOR, ZRSR2, IDH2, SETBP1, DDX41, ETV6, ETNK1

mutations <- c("ASXL1", "SRSF2", "DNMT3A", "RUNX1", "U2AF1", "EZH2", 
    "BCOR", "ZRSR2", "IDH2", "SETBP1", "DDX41", "ETV6", "ETNK1")

png("figures/GESMD_IWS_clustering/subgroup_exploration/ASXL1_subgroups.png", width = 2000, height = 1000, res = 300)
rbind(gesmd %>% select(sub_group, ASXL1) %>% mutate(dataset = "GESMD"),
      clinical_blasts %>% select(sub_group, ASXL1) %>% mutate(dataset = "IWS")) %>%
      change_levels() %>%
      group_by(sub_group, dataset) %>%
      summarize(Freq = mean(ASXL1)*100) %>%
      ggplot(aes(x = dataset, y = Freq, fill = sub_group)) +
      geom_bar(stat = "identity") +
      facet_grid(~ sub_group) +
        scale_fill +
      ylab("Frequency (%)") +
      xlab("") +
      ggtitle("ASXL1 frequency") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_blank(), strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5))
dev.off()

plot_gene_frequency <- function(gene, scale_fill, output_dir = "figures/GESMD_IWS_clustering/subgroup_exploration") {
    png_filename <- file.path(output_dir, paste0(gene, "_subgroups.png"))
    rbind(
        gesmd %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "GESMD"),
        clinical_blasts %>% select(sub_group, !!sym(gene)) %>% mutate(dataset = "IWS")
    ) %>%
        change_levels() %>%
        group_by(sub_group, dataset) %>%
        summarize(Freq = mean(!!sym(gene)) * 100, .groups = "drop") %>%
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
lapply(mutations[mutations != "DDX41"], plot_gene_frequency, scale_fill = scale_fill)

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
            N = n(),
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
            ASXL1 = getProp(ASXL1 == 1, ASXL1),
            SRSF2 = getProp(SRSF2 == 1, SRSF2),
            DNMT3A = getProp(DNMT3A == 1, DNMT3A),
            RUNX1 = getProp(RUNX1 == 1, RUNX1),
            U2AF1 = getProp(U2AF1 == 1, U2AF1),
            EZH2 = getProp(EZH2 == 1, EZH2),
            BCOR = getProp(BCOR == 1, BCOR),
            ZRSR2 = getProp(ZRSR2 == 1, ZRSR2),
            IDH2 = getProp(IDH2 == 1, IDH2),
            SETBP1 = getProp(SETBP1 == 1, SETBP1),
            ETV6 = getProp(ETV6 == 1, ETV6),
            ETNK1 = getProp(ETNK1 == 1, ETNK1),
            DDX41 = getProp(DDX41 == 1, DDX41)    
           )
}

gesmd_mod <- gesmd %>%
    mutate(IPSSM = ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_ALTO", "High",
        ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_BAJO", "Low", 
           ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_ALTO", "Moderate-High",
            ifelse(IPSSM == "IPSS_MOL_GRUPO_RIESGO_MOD_BAJO", "Moderate-Low", IPSSM)) )),
        DDX41 = NA) %>%
    change_levels()
gesmd_subgroup <- rbind(gesmd_mod  %>%
    group_by(sub_group) %>%
    summarize_fun(),
    cbind(tibble(sub_group = "Total"),  summarize_fun(gesmd_mod))
) %>%
    t() 
gesmd_subgroup_df <- as.data.frame(gesmd_subgroup[-1, ]) 
colnames(gesmd_subgroup_df) <- gesmd_subgroup[1, ]
rownames(gesmd_subgroup_df) <- rownames(gesmd_subgroup)[-1]


## Tests
poisson_test <- function(var, df){
    poiss_base <- glm(formula (paste(var, " ~ 1")), df, 
                      family = "poisson")
    poiss_lm <- glm(formula (paste(var, " ~ sub_group")), df, 
                      family = "poisson")
    anova_res <- anova(poiss_base, poiss_lm, test = "Chisq")                      
    anova_res$`Pr(>Chi)`[2]
}

lm_test <- function(var, df){
    lm_res <- summary(lm(formula (paste(var, " ~ sub_group")), df))
    fstats <- lm_res$fstatistic
    pf(fstats[1], fstats[2], fstats[3], lower.tail = FALSE)
 }

chisq_test <- function(var, df){
    chisq_res <- chisq.test(table(df[[var]], df$sub_group))
    chisq_res$p.value
}

# Perform tests
gesmd_test <- c(sapply(c("SEX", "consensus", "IPSSM", mutations[mutations != "DDX41"]), chisq_test, df = gesmd_mod),
            sapply(c("AGE", "HB", "PLT"), lm_test, df = gesmd_mod),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = gesmd_mod))
names(gesmd_test) <- c("Females", "Low blasts", "IPSSM_NA", mutations[mutations != "DDX41"], "Age", "HB", "PLT", 
    "BM Blasts", "WBC count", "Neutrophil Count", "Monocyte Count") 

gesmd_subgroup_df$p_value <- gesmd_test[rownames(gesmd_subgroup_df)]

write.table(gesmd_subgroup_df, 
            file = "results/GESMD_IWS_clustering/gesmd_subgroup_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)

### IWS
iws_subgroup <- rbind(clinical_blasts  %>%
    group_by(sub_group) %>%
    summarize_fun(),
    cbind(tibble(sub_group = "Total"),  summarize_fun(clinical_blasts))
) %>%
    t() 
iws_subgroup_df <- as.data.frame(iws_subgroup[-1, ]) 
colnames(iws_subgroup_df) <- iws_subgroup[1, ]
rownames(iws_subgroup_df) <- rownames(iws_subgroup)[-1]



# Perform tests
iws_test <- c(sapply(c("SEX", "consensus", "IPSSM", mutations), chisq_test, df = clinical_blasts),
            sapply(c("AGE", "HB", "PLT"), lm_test, df = clinical_blasts),
            sapply(c("BM_BLAST", "WBC", "ANC", "MONOCYTES"), poisson_test, df = clinical_blasts))
names(iws_test) <- c("Females", "Low blasts", "IPSSM_NA", mutations, "Age", "HB", "PLT", 
    "BM Blasts", "WBC count", "Neutrophil Count", "Monocyte Count") 

iws_subgroup_df$p_value <- iws_test[rownames(iws_subgroup_df)]

write.table(iws_subgroup_df, 
            file = "results/GESMD_IWS_clustering/IWS_subgroup_descriptives.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)
