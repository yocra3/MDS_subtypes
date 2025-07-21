library(tidyverse)
library(stringr)
library(ggpubr)

scores <- read_table("results/gnn/full_input_model/v1_boolean/v1_boolean_individual_gene_predictions.csv") %>%
  select(-Combination)
scores_per_gene <- read_table("results/gnn/full_input_model/v1_boolean/v1_boolean_individual_gene_predictions_per_gene.csv")

base <- subset(scores, Gene == "base")
hist(base$Score)


scores_mode <- scores %>%
  group_by(ID) %>%
  mutate(score_diff = Score - Score[Gene == "base"]) %>%
  ungroup() %>%
  mutate(n_genes = str_count(Gene, "_") + 1)
  
single_genes %>% 
  group_by(Gene) %>% 
  summarize(m = median(score_diff)) %>%
  arrange(desc(m))

single_genes %>% 
  group_by(Gene) %>% 
  summarize(m = median(score_diff)) %>%
  arrange(m)


single_genes <- subset(scores_mode, n_genes == 1 & Gene != "base")
single_genes %>%
  group_by(Gene) %>%
  summarize(m = median(score_diff)) %>%
  arrange(desc(m))

boxplot(single_genes$score_diff ~ single_genes$Gene, las = 2, main = "Single Gene Predictions", ylab = "Score Difference from Base")

tapply(single_genes$score_diff, single_genes$Gene,  summary)

get_combo_name <- function(gvec) {
  if (length(gvec) == 0) return("base")
  paste(sort(gvec), collapse = "_")
}

get_scores <- function(df){
  f_values <- setNames(df$Score, df$Gene)
  names(f_values) <- sapply(names(f_values), function(x){
    vec <- strsplit(x, "_")[[1]]
    get_combo_name(vec)
  })
  f_values
}
fac <- factorial

compute_Shap_vals <- function(df_ind){
  
  
  scores <- get_scores(df_ind)
  
  genes <- df_ind %>%
    filter(n_genes == 1, Gene != "base") %>%
    pull(Gene) %>%
    unique()
  
  n <- length(genes)
  shap_vals <- sapply(genes, function(i){
    
    others <- setdiff(genes, i)
    
    # 1) Caso S = ∅
    weight0 <- fac(0) * fac(n - 1) / fac(n)  # = (n-1)! / n! = 1/n
    f_empty <- scores["base"]
    f_i     <- scores[i]
    shap_val <- weight0 * (f_i - f_empty)
    
    # 2) Casos S ⊂ others, |S| ≥ 1
    for (s in seq_along(others)) {
      combs <- combn(others, s, simplify = FALSE)
      for (S in combs) {
        S_name  <- get_combo_name(S)
        Si_name <- get_combo_name(c(S, i))
        
        fS  <- scores[S_name]
        fSi <- scores[Si_name]
        
        weight <- fac(length(S)) * fac(n - length(S) - 1) / fac(n)
        shap_val <- shap_val + weight * (fSi - fS)
      }
    }
    
    shap_val
  })
  names(shap_vals) <- genes
  shap_vals
  }

patients <- unique(scores_mode$ID)
shaps <- lapply(patients, function(id){
  df_ind <- subset(scores_mode, ID == id)
  vals <- compute_Shap_vals(df_ind)
  tibble(ID = id, Gene = names(vals), SHAP = vals)
}) %>%
  Reduce(f = rbind)

IPSSM_genes <- c("TP53", "MLL", "FLT3", "SF3B1", "NPM1", "RUNX1", "NRAS", "ETV6",
                 "IDH2", "CBL", "EZH2", "U2AF1", "SRSF2", "DNMT3A", "ASXL1", 
                 "KRAS", "BCOR", "BCORL1", "CEBPA", "ETNK1", "GATA2", "GNB1", 
                 "IDH1", "NF1", "PHF6", "PPM1D", "PRPF8", "PTPN11", "SETBP1", 
                 "STAG2", "WT1")

gene_tab <- shaps %>% group_by(Gene) %>%
  summarize(n = n()) %>%
  mutate(IPSSM = ifelse(Gene %in% IPSSM_genes, "IPSSM", "Other"),
         Frequency = ifelse(n < 10, "10-", "10+"))

gene_order <-  shaps %>% 
  group_by(Gene) %>% 
  summarize(m = median(SHAP)) %>% 
  arrange(desc(m)) %>% 
  pull(Gene)

png("figures/gnn/SHAP_genes.png", width = 800, height = 600)
shaps %>%
  left_join(gene_tab, by = "Gene")  %>%
  mutate(Gene = factor(Gene, levels = gene_order)) %>%
  ggplot(aes(x = Gene, y = SHAP, fill = IPSSM)) +
  geom_boxplot() +
  labs(title = "SHAP Values by Gene", y = "SHAP Value", x = "Gene") +
  theme_bw() +
  facet_wrap(Frequency ~ ., scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

## Compute shap values for models per gene
tested_genes <- unique(scores_per_gene$Testing_Gene)
scores_per_gene_mode <- scores_per_gene %>%
  mutate(Gene = Gene_Combination) %>%
  group_by(ID) %>%
  mutate(score_diff = Score - Score[Gene == "base"]) %>%
  ungroup() %>%
  mutate(n_genes = str_count(Gene, "_") + 1)
  

shaps_per_gene <- lapply(tested_genes, function(gene){
  df_gene <- subset(scores_per_gene_mode, Testing_Gene == gene)
  patients <- unique(df_gene$ID)
  out_df <- lapply(patients, function(id){
    df_ind <- subset(df_gene, ID == id)
    vals <- compute_Shap_vals(df_ind)
    tibble(ID = id, Gene = names(vals), SHAP = vals)
  }) %>%
    Reduce(f = rbind) %>%
    mutate(Testing_Gene = gene)
}) %>%
  Reduce(f = rbind)

comb_shap <- filter(shaps_per_gene, Gene == Testing_Gene) %>%
  left_join(shaps, by = c("ID", "Gene"), suffix = c("_gene", "_full")) %>%
  mutate(SHAP_diff = SHAP_gene - SHAP_full)

png("figures/gnn/SHAP_comparison.png", width = 800, height = 600)
ggplot(comb_shap, aes(x = Gene, y = SHAP_diff)) +
  geom_boxplot() +
  labs(title = "SHAP Value Differences by Gene", y = "SHAP Difference", x = "Gene") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()



## Explore factors affecting SHAP values
reduced_embeddings <- read_delim("results/preprocess/gene_embedding_reduced.tsv", delim = "\t")

emb_mat <- as.matrix(reduced_embeddings[, -1])
rownames(emb_mat) <- unlist(reduced_embeddings[, 1])

filt_mat <- emb_mat[rownames(emb_mat) %in% shaps_per_gene$Gene, ]
emb_cor <- cor(t(filt_mat), use = "pairwise.complete.obs")
diag(emb_cor) <- 0  # Remove self-correlation

gene_cor <- data.frame(Gene = rownames(emb_cor), 
                   cor = matrixStats::colMaxs(emb_cor, na.rm = TRUE),
                   redundancy = apply(emb_cor, 1, function(x) sum(x > 0.7, na.rm = TRUE)))


patient_per_gene <- comb_shap %>%
  group_by(Gene) %>%
  summarize(n_patients = n_distinct(ID)) %>%
  ungroup()

c_index <- read_delim("results/gnn/full_input_model_per_gene/v1_boolean_gene_summary_all.csv")

comb_shap_sum <- comb_shap %>%
  group_by(Testing_Gene) %>%
  summarize(median_diff = median(SHAP_diff, na.rm = TRUE),
        SHAP_cor = cor(SHAP_gene, SHAP_full, use = "complete.obs")) %>%
  left_join(gene_cor, by = c("Testing_Gene" = "Gene")) %>%
  left_join(patient_per_gene, by = c("Testing_Gene" = "Gene")) %>%
  left_join(c_index %>% select(gene, best_val_c_index), by = c("Testing_Gene" = "gene")) %>%
  mutate(IPSSM = ifelse(Testing_Gene %in% IPSSM_genes, "IPSSM", "Other")) %>%
  arrange(desc(SHAP_cor)) 

summary(lm(SHAP_cor ~ cor + redundancy + redundancy + best_val_c_index, comb_shap_sum))


png("figures/gnn/SHAP_genes_gene_model.png", width = 800, height = 600)
filter(shaps_per_gene, Gene == Testing_Gene) %>%
  left_join(gene_tab, by = "Gene")  %>%
  mutate(Gene = factor(Gene, levels = gene_order)) %>%
  ggplot(aes(x = Gene, y = SHAP, fill = IPSSM)) +
  geom_boxplot() +
  labs(title = "SHAP Values by Gene (gene model)", y = "SHAP Value", x = "Gene") +
  theme_bw() +
  facet_wrap(Frequency ~ ., scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()


png("figures/gnn/SHAP_comparison_scatter.png", width = 2000, height = 1600)
comb_shap %>%
  mutate(Testing_Gene = factor(Testing_Gene, levels = comb_shap_sum$Testing_Gene)) %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full)) +
  geom_point() +
  labs(title = "SHAP Value Differences by Gene", y = "SHAP Full", x = "SHAP Gene") +
  theme_bw() +
  facet_wrap(~ Testing_Gene, scales = "free") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")
dev.off()

## Define scores for gene combinations

scores_mode2 <- scores_mode %>%
  filter(n_genes > 1) %>%
  rowwise() %>%
  mutate(
    genes_split = str_split(Gene, "_")[1],
    score_diff_ind_genes = sum(single_genes$score_diff[
      single_genes$ID == ID & single_genes$Gene %in% genes_split
    ], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-genes_split) %>%
  mutate(inter_effect = score_diff - score_diff_ind_genes)

png("figures/gnn/SHAP_inter_effect.png", width = 1500, height = 1100, res = 300)
ggplot(scores_mode2, aes(x = factor(n_genes), y = inter_effect)) +
  geom_boxplot() +
  labs(title = "SHAP Interaction Effect by Number of Genes", y = "Interaction Effect", x = "Number of Genes") +
  theme_bw()
dev.off()


shap_5q <- mutate(shaps, del5q = ID %in% subset(clinical, del5q == 1)$ID)

inters <- filter(scores_mode2, abs(inter_effect) > log2(1.5))
sort(table(unlist(strsplit(inters$Gene, "_"))))

pos_int_ids <- filter(scores_mode2, inter_effect > 0.1) %>% pull(ID)
scores_mode2 %>%
  filter(ID %in% pos_int_ids) %>%
  group_by(ID) %>%
  summarize(median_inter_effect = median(inter_effect, na.rm = TRUE),
  mean_pos = mean(inter_effect > 0)) %>%
  arrange(desc(median_inter_effect)) 

sel_pos <- scores_mode2 %>%
  filter(ID %in% pos_int_ids) %>%
  group_by(ID) %>%
  summarize(median_inter_effect = median(inter_effect, na.rm = TRUE),
  mean_pos = mean(inter_effect > 0)) %>%
  arrange(desc(median_inter_effect)) %>%
  filter(mean_pos < 0.5) %>%
  pull(ID)

filter(scores_mode2, ID == "E-H-122335-T1-1-D1-1") %>%
  data.frame()
filter(single_genes, ID == "E-H-122335-T1-1-D1-1") %>%
  data.frame()

filter(scores_mode2, ID == "E-H-110804-T1-1-D1-1") %>%
  data.frame()
filter(single_genes, ID == "E-H-110804-T1-1-D1-1") %>%
  data.frame()


filter(scores_mode2, ID == "E-H-118502-T1-1-D1-1") %>%
  data.frame()
filter(single_genes, ID == "E-H-118502-T1-1-D1-1") %>%
  data.frame()

filter(scores_mode2, ID == "E-H-100116-T1-1-D1-1") %>%
  data.frame()
filter(single_genes, ID == "E-H-100116-T1-1-D1-1") %>%
  data.frame()


filter(scores_mode2, ID == "E-H-110780-T1-1-D1-1") %>%
  data.frame()

cot <- scores_mode2 %>%
  filter(ID %in% pos_int_ids) %>%
  group_by(ID) %>%
  mutate(SF3B1 = ifelse(str_detect(Gene, "SF3B1"), 1, 0)) %>%
   summarize(mean_pos_sf3b1 = sum(SF3B1 & inter_effect > 0)/sum(SF3B1),
             mean_sf3b1_pos = sum(SF3B1 & inter_effect > 0)/sum(inter_effect > 0)) 

sum(cot$mean_sf3b1_pos > 0)
sum(filter(single_genes, Gene == "SF3B1")$ID %in% scores_mode2$ID)