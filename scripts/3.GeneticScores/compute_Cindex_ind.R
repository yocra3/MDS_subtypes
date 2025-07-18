#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compute individual C-index for IPSSM and the GNN models
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
load("results/preprocess/clinical_preproc.Rdata")

clin_vals <- read_delim("results/gnn/all_patients_predictions.tsv", delim = "\t")

clin_all <- left_join(clin_vals %>% select(gnn_bool_optim_pred, gnn_vaf_optim_pred, ID),
    clinical, by = "ID")

compute_ind_Cindex <- function(df, score, time, event){

    samps <- nrow(df)
    risk_row <- matrix(df[[score]], nrow = samps, ncol = samps)
    risk_col <- matrix(df[[score]], nrow = samps, ncol = samps, byrow = TRUE)
    risk_bool <- risk_row > risk_col
    risk_bool2 <- risk_row == risk_col

    time_row <- matrix(df[[time]], nrow = samps, ncol = samps)
    time_col <- matrix(df[[time]], nrow = samps, ncol = samps, byrow = TRUE)
    time_bool <- time_row < time_col

    event_row <- matrix(df[[event]], nrow = samps, ncol = samps)
    event_col <- matrix(df[[event]], nrow = samps, ncol = samps, byrow = TRUE)
    valid <- time_row < time_col & event_row == 1 | time_row > time_col & event_col == 1   

    c_index_mat <- (risk_bool == time_bool)
    c_index_mat[(risk_bool2 * valid) == 1] <- 0.5
    c_index_mat[!valid] <- NA
    rownames(c_index_mat) <- df$ID

    rowMeans(c_index_mat, na.rm = TRUE)
}
clin_all$gnn_boolean_cind <- compute_ind_Cindex(clin_all, "gnn_bool_optim_pred", "LFS_YEARS", "LFS_STATUS")
clin_all$gnn_vaf_cind <- compute_ind_Cindex(clin_all, "gnn_vaf_optim_pred", "LFS_YEARS", "LFS_STATUS")
clin_all$IPSSM_cind <- compute_ind_Cindex(clin_all, "IPSSM_SCORE", "LFS_YEARS", "LFS_STATUS")
clin_all$gnn_boolean_cind_diff <- clin_all$gnn_boolean_cind - clin_all$IPSSM_cind
clin_all$gnn_boolean_cut <- cut(clin_all$gnn_boolean_cind_diff, breaks = seq(-0.6, 0.7, 0.1))
clin_all$gnn_vaf_cind_diff <- clin_all$gnn_vaf_cind - clin_all$IPSSM_cind
clin_all$gnn_vaf_cut <- cut(clin_all$gnn_vaf_cind_diff, breaks = seq(-0.6, 0.7, 0.1))

clin_all$gnn_boolean_z <- scale(clin_all$gnn_bool_optim_pred )
clin_all$gnn_vaf_z <- scale(clin_all$gnn_vaf_optim_pred )

cor(clin_all$gnn_bool_optim_pred, clin_all$IPSSM_SCORE)
cor(clin_all$gnn_boolean_cind, clin_all$IPSSM_cind, use = "complete.obs")

png("figures/gnn_models/comp_boolean_ipssm.png", res = 300, width = 1500, height = 1000)
ggplot(clin_all, aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE)) +
    geom_point(alpha = 0.1) + 
    theme_bw()
dev.off()

patient_feats <- c("SEX", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex")
names(patient_feats) <- patient_feats
lapply(patient_feats, function(feat){
    fisher.test(table(clin_all$gnn_boolean_cind_diff < -0.1, clin_all[[feat]]))
})
table(clin_all$gnn_boolean_cind_diff < -0.1, clin_all$SEX)
table(clin_all$gnn_boolean_cind_diff < -0.1, clin_all$delY)

lapply(patient_feats, function(feat){
    p <- clin_all %>%  
        arrange(.data[[feat]]) %>%  
        ggplot(aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE, color = factor(.data[[feat]] ))) +
        geom_point() + 
        scale_color_manual(values = c("black", "#05d1d1")) + 
        ggtitle(feat) +
        xlab("C-index difference") +
        ylab("Score difference") +
        theme_bw()
    ggsave(paste0("figures/gnn_models/comp_boolean_ipssm_", feat, ".png"), p, units = "px", 
    width = 1500, height = 1000)
})
lapply(patient_feats, function(feat){
    summary(lm(formula(paste("gnn_boolean_cind_diff ~ ", feat)), clin_all))
})
genes <- colnames(clin_all)[323:403]

lapply(genes, function(gene){
    p <- clin_all %>%  
        arrange(.data[[gene]]) %>%  
        ggplot(aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE, color = factor(.data[[gene]] ))) +
        geom_point() + 
        scale_color_manual(values = c("black", "#05d1d1")) + 
        ggtitle(gene) +
        xlab("C-index difference") +
        ylab("Score difference") +
        theme_bw()
    ggsave(paste0("figures/gnn_models/comp_boolean_ipssm_", gene, ".png"), p, units = "px", 
    width = 1500, height = 1000)
})
genes_lm <- lapply(sort(genes), function(feat){
    summary(lm(formula(paste("gnn_boolean_cind_diff ~ ", feat)), clin_all))
})
names(genes_lm) <- sort(genes)
genes_p <- sapply(genes_lm, function(x) x$coefficients[2, 4])

genes_lm[names(genes_p)[genes_p < 0.05]]
## Conclusiones comparación 1
#' - DelY, DDX41, GNAS: más riesgo en modelo GNN. A veces es mejor predictor, otras no.
#' - Complex: en general, peor predicción para GNN. Relación dificil
#' - TP53: no hay correlación clara. Los multi tienen más riesgo en IPSSM y los mono en GNN. La predicción es parecida
#' - del5q, ASXL1, CBL, ARID2, ASXL2, ATRX, BRAF, BRCC3, CREBBP, DDX23, EP300, FLT3, IRF1, KMT2C, PHF6, PPM1D, PRPF8, PTPN11, RB1, ROBO1, WT1: Predicción mejor en GNN
#' - del7, del7q, del20q, 8+, DNMT3A, KRAS, NPM1, NRAS, RUNX1, SRSF2, BCOR, CALR, CSF3R, CTCF, CUX1, DDX54, EED, ETNK1, GATA2, GNB1, IDH1, KIT; KMT2D, MGA, NF1, NPM1, PIK3CA, RAD50, ROBO2, SH2B3, SMC3, STAG2, STAT3, TERT, TET2 (bi y mono), U2AF2, ZRSR2: pocas diferencias
#' - ETV6, MLL_PTD: más riesgo IPSSM. Predicción parecida
#' - EZH2, IDH2, U2AF1, BCORL1, CEBPA, CSKNK1A, DDX4, FLT3_ITD, JAK2, KDM6A, LUC7L2, MPL, SETBP1, SETD2: predicción mejor en IPSSM
#' - NFE2, RAD21, SF3B1, SMC1A, SUZ12, ZBTB33: más riesgo GNN y mejor predicción IPSSM
#' KDM6A: peor predicción GNN

png("figures/gnn_models/comp_boolean_ipssm_ngenes.png", res = 300, width = 1500, height = 1000)
     clin_all %>%  
        ggplot(aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE, color = factor(N_mutations))) +
        geom_point() + 
        ggtitle("N mutations") +
        xlab("C-index difference") +
        ylab("Score difference") +
        theme_bw()
dev.off()

## No hay correlación entre el número de mutaciones y el performance
summary(lm(gnn_boolean_cind_diff ~ N_mutations, clin_all))

clin_all$SF3B1_del5q <- clin_all$SF3B1 * clin_all$del5q
clin_all$SF3B1_delY <- clin_all$SF3B1 * clin_all$delY
clin_all$SETBP1_del7 <- clin_all$SETBP1 * clin_all$del7

summary(lm(gnn_boolean_cind_diff ~ SF3B1_delY, clin_all))
summary(lm(gnn_boolean_cind_diff ~ SF3B1_del5q, clin_all))
summary(lm(gnn_boolean_cind_diff ~ SETBP1_del7, clin_all))


summary(lm(gnn_boolean_cind_diff ~ CYTO_IPSSR, clin_all))
png("figures/gnn_models/comp_boolean_ipssm_CYTO_IPSSR.png", res = 300, width = 1500, height = 1000)
     clin_all %>%  
        ggplot(aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE, color = factor(CYTO_IPSSR))) +
        geom_point() + 
        ggtitle("CYTO_IPSSR") +
        xlab("C-index difference") +
        ylab("Score difference") +
        theme_bw()
dev.off()

png("figures/gnn_models/comp_boolean_ipssm_del7_SETBP1.png", res = 300, width = 1500, height = 1000)
     clin_all %>%  
        arrange(SETBP1_del7) %>%
        ggplot(aes(x = gnn_boolean_cind - IPSSM_cind, y = gnn_boolean_z - IPSSM_SCORE, color = factor(SETBP1_del7))) +
        geom_point() + 
        scale_color_manual(values = c("black", "#05d1d1")) + 
        ggtitle("SETBP1 del7") +
        xlab("C-index difference") +
        ylab("Score difference") +
        theme_bw()
dev.off()

png("figures/gnn_models/ind_cindex_boolean_ipssm.png", res = 300, width = 1500, height = 1000)
ggplot(clin_vals, aes(x = IPSSM_cind, y = gnn_boolean_cind, color = gnn_bool_optim_pred)) +
    geom_point() + 
    scale_color_viridis_c() +
    theme_bw()
dev.off()


png("figures/gnn_models/ind_cindex_vaf_ipssm.png", res = 300, width = 1500, height = 1000)
ggplot(clin_vals, aes(x = IPSSM_cind, y = gnn_vaf_cind, color = gnn_vaf_optim_pred)) +
    geom_point() + 
    scale_color_viridis_c() +
    theme_bw()
dev.off()

## Select patients most failing in GNN boolean with respect to IPSSM
fail_gnn_bool <- subset(clin_vals, gnn_boolean_cind - IPSSM_cind < -0.1)
fail_gnn_vaf <- subset(clin_vals, gnn_vaf_cind - IPSSM_cind < -0.1)

overlaps <- function(x, y){
    lengths(list(common = intersect(x, y), bool = setdiff(x, y), vaf = setdiff(y, x)))
}
overlaps(fail_gnn_bool$ID, fail_gnn_vaf$ID)


fail2_gnn_bool <- subset(clin_vals, gnn_boolean_cind - IPSSM_cind < -0.3)

## Load mutation data and select relevant columns


max_n <-  maf %>% group_by(ID, GENE) %>% 
    summarize(n = n()) %>% group_by(ID) %>% 
    summarize(max_n = max(n)) %>%
    mutate(bad_pat_bool = ifelse(ID %in% fail_gnn_bool$ID, "Bad", "Good"), 
    bad_pat_vaf = ifelse(ID %in% fail_gnn_vaf$ID, "Bad", "Good"))

table(max_n$max_n, max_n$bad_pat_bool)
table(max_n$max_n, max_n$bad_pat_vaf)

n_genes <- maf %>% group_by(ID) %>% 
    summarize(n_genes = length(unique(GENE)))  %>%
    mutate(bad_pat_bool = ifelse(ID %in% fail_gnn_bool$ID, "Bad", "Good"), 
    bad_pat_vaf = ifelse(ID %in% fail_gnn_vaf$ID, "Bad", "Good"),
    n_genes_cat = ifelse(n_genes > 4, "5+", n_genes))


table(n_genes$n_genes, n_genes$bad_pat_bool)
table(n_genes$n_genes, n_genes$bad_pat_vaf)

chisq.test(table(n_genes$n_genes, n_genes$bad_pat_bool))
chisq.test(table(n_genes$n_genes, n_genes$bad_pat_vaf))

table(n_genes$n_genes_cat, n_genes$bad_pat_bool)
table(n_genes$n_genes_cat, n_genes$bad_pat_vaf)

chisq.test(table(n_genes$n_genes_cat, n_genes$bad_pat_bool))
chisq.test(table(n_genes$n_genes_cat, n_genes$bad_pat_vaf))


png("figures/gnn_models/gnn_boolean_badpat_vaf.png", res = 300, width = 5000, height = 4000)
 maf %>% filter(GENE %in% 
    (group_by(maf, GENE) %>% summarize(bad = any(bad_pat_bool == "Bad")) %>% filter(bad) %>% pull(GENE))
 ) %>%
ggplot(aes(x = bad_pat_bool, y = VAF)) +
    geom_boxplot() + 
    facet_wrap(~ GENE) +
    theme_bw()
dev.off()

bool_gene_bad_sum <- maf %>% filter(GENE %in% 
    (group_by(maf, GENE) %>% summarize(bad = any(bad_pat_bool == "Bad")) %>% filter(bad) %>% pull(GENE))
 ) %>%
    group_by(bad_pat_bool, GENE) %>%
    summarize(n = n()) %>% 
    pivot_wider(names_from = bad_pat_bool, values_from = n) %>%
    mutate(bad_prop = Bad/(Bad + Good)) %>%
    arrange(desc(bad_prop)) 

bad_genes_bool <- filter(bool_gene_bad_sum, bad_prop > nrow(fail_gnn_bool) / nrow(clin_vals) & Bad > 5)$GENE

maf_bool_wide <- maf %>% select(ID, GENE, bad_pat_bool) %>% 
    distinct() %>%
    mutate(status = 1) %>%
    pivot_wider(names_from = GENE, values_from = status, values_fill = 0)


maf_bool_wide_comb <- maf_bool_wide %>%
  rowwise() %>%
  mutate(
    gene_comb = {
      vars <- names(cur_data())[-c(1:2)]
      paste(vars[c_across(3:last_col()) == 1], collapse = "; ")
    }
  ) %>%
  ungroup()

bool_gene_bad_sum$enrich_pval <- sapply(bool_gene_bad_sum$GENE, function(gene){
    fisher.test(table(maf_bool_wide[[gene]], maf_bool_wide$bad_pat_bool))$p.value
})

bool_gene_bad_sum$vaf_pval <- sapply(bool_gene_bad_sum$GENE, function(gene){
    if (gene %in% c("MLL", "PRPF40A")){
        return(NA)
    } else {
        summary(lm(VAF ~ bad_pat_bool, data = maf, subset = GENE == gene))$coefficients[2, 4]
    }
})

filter(bool_gene_bad_sum, (enrich_pval < 0.05 | vaf_pval < 0.05) & Bad > 5)


png("figures/gnn_models/gnn_vaf_badpat_vaf.png", res = 300, width = 5000, height = 4000)
 maf %>% filter(GENE %in% 
    (group_by(maf, GENE) %>% summarize(bad = any(bad_pat_vaf == "Bad")) %>% filter(bad) %>% pull(GENE))
 ) %>%
ggplot(aes(x = bad_pat_vaf, y = VAF)) +
    geom_boxplot() + 
    facet_wrap(~ GENE) +
    theme_bw()
dev.off()

vaf_gene_bad_sum <- maf %>% filter(GENE %in% 
    (group_by(maf, GENE) %>% summarize(bad = any(bad_pat_vaf == "Bad")) %>% filter(bad) %>% pull(GENE))
 ) %>%
    group_by(bad_pat_vaf, GENE) %>%
    summarize(n = n()) %>% 
    pivot_wider(names_from = bad_pat_vaf, values_from = n) %>%
    mutate(bad_prop = Bad/(Bad + Good)) %>%
    arrange(desc(bad_prop)) 


maf_vaf_wide <- maf %>% select(ID, GENE, bad_pat_vaf) %>% 
    distinct() %>%
    mutate(status = 1) %>%
    pivot_wider(names_from = GENE, values_from = status, values_fill = 0)

vaf_gene_bad_sum$enrich_pval <- sapply(vaf_gene_bad_sum$GENE, function(gene){
    fisher.test(table(maf_vaf_wide[[gene]], maf_vaf_wide$bad_pat_vaf))$p.value
})

vaf_gene_bad_sum$vaf_pval <- sapply(vaf_gene_bad_sum$GENE, function(gene){
    if (gene %in% c("MLL", "PRPF40A")){
        return(NA)
    } else {
        summary(lm(VAF ~ bad_pat_vaf, data = maf, subset = GENE == gene))$coefficients[2, 4]
    }
})

filter(vaf_gene_bad_sum, (enrich_pval < 0.05 | vaf_pval < 0.05) & Bad > 5)


gene_tab <- maf %>% group_by(ID, bad_pat_bool) %>% summarize(n = sum(GENE == "EZH2"))
chisq.test(table(gene_tab$bad_pat_bool, gene_tab$n))
chisq.test(table(gene_tab$bad_pat_bool, gene_tab$n)[, -1])

## Higher number of mutations for SF3B1