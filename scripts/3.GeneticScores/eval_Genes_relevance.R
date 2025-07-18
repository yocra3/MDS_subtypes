library(tidyverse)
library(stringr)

scores <- read_table("results/v1_boolean_individual_gene_predictions.csv") %>%
  select(-Combination)

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

shaps %>%
  left_join(gene_tab, by = "Gene")  %>%
  ggplot(aes(x = Gene, y = SHAP, fill = IPSSM)) +
  geom_boxplot() +
  labs(title = "SHAP Values by Gene", y = "SHAP Value", x = "Gene") +
  theme_bw() +
  facet_wrap(Frequency ~ ., scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(angle = 90))

library(tabulapdf)

# Extrae todas las tablas de la página 1
tables <- extract_tables("table4_evidoa2200008_appendix.pdf")


## Define scores fr
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
hist(scores_mode2$inter_effect)
boxplot(inter_effect ~ n_genes, scores_mode2)

table(ifelse(scores_mode2$inter_effect > 0.585, "positive",
             ifelse(scores_mode2$inter_effect < -0.585, "negative", "neutral")),
      scores_mode2$n_genes)

phi_int <- matrix(0, nrow=n, ncol=n, dimnames=list(genes, genes))

for(i in genes){
  for(j in genes){
    if(j<=i) next
    others <- setdiff(genes, c(i,j))
    for(s in seq_along(others)){
      for(S in combn(others, s, simplify=FALSE)){
        Sname <- get_combo_name(S)
        fS <- f_values[Sname]
        fSi <- f_values[get_combo_name(c(S,i))]
        fSj <- f_values[get_combo_name(c(S,j))]
        fSij<- f_values[get_combo_name(c(S,i,j))]
        w <- fac(s)*fac(n-s-2)/(2*fac(n-1))
        phi_int[i,j] <- phi_int[i,j] + w*(fSij - fSi - fSj + fS)
      }
    }
  }
}


shap_5q <- mutate(shaps, del5q = ID %in% subset(clinical, del5q == 1)$ID)

inters <- filter(scores_mode2, abs(inter_effect) > log2(1.5))
sort(table(unlist(strsplit(inters$Gene, "_"))))
