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

## Define normalized scores
final_scores <- scores_mode %>% 
  group_by(ID) %>% 
  filter(n_genes == max(n_genes)) %>%
  mutate(nrow = n()) %>%
  filter(!(Gene == "base" & nrow == 2)) %>%
  select(ID, Score)
write.table(final_scores, "results/gnn/full_input_model/v1_boolean/v1_boolean_normalized_scores.tsv", 
  sep = "\t", row.names = FALSE, quote = FALSE)


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

## Explore SHAP values in relevant genes
shap_patient <- left_join(shaps,
  clinical %>% select(ID, any_of(IPSSM_genes), TP53multi, TP53mono, 
    TET2bi, TET2other, N_mutations, complex, del5q, plus8, del20q, del7, del7q, delY
    CYTO_IPSSR), by = "ID")
maf <- read_delim("data/IPSSMol/maf.tsv")

isMain <- function(gene, pat, maf) {
  patient <- subset(maf, ID == pat)
  max_vaf <- max(patient$VAF, na.rm = TRUE)
  gene_maf <- max(subset(patient, GENE == gene)$VAF)
  return (gene_maf + 0.1 > max_vaf)
}

### TET2
tet2_f <- shap_patient %>%
  filter(Gene == "TET2") %>%
  mutate(HR = 2**SHAP,
    TET2_Status = ifelse(TET2bi == 1, "Bi-allelic", "Mono-allelic"),
    Context = ifelse(SF3B1 == 1, "TET2 + SF3B1", 
      ifelse(RUNX1 == 1, "TET2 + RUNX1", 
        ifelse(CBL == 1, "TET2 + CBL", 
        ifelse(SRSF2 == 1, "TET2 + SRSF2", 
            ifelse(complex == 1, "TET2 + complex", 
              ifelse(del7 == 1 | del7q == 1, "TET2 + del7/del7q", 
                    ifelse(N_mutations == 1, "Only TET2", 
                      ifelse(N_mutations == 2, "TET2 and 1 mut.", "TET2 and 2+ mut."))))))))) %>%
    rowwise() %>%
    mutate(main_mut = isMain("U2AF1", ID, maf)) 
summary(lm(SHAP ~ ., select(tet2_f, -c(ID, TET2bi, Gene, Context, HR))))
summary(lm(SHAP ~ SF3B1 + RUNX1 + CBL + SRSF2 + STAG2 + TET2other + complex  + del7 + del7q + N_mutations + main_mut, tet2_f))

png("figures/gnn/SHAP_full_TET2.png", width = 2500, height = 1500, res = 300)
ggplot(tet2_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for TET2", y = "HR", x = "Context") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ TET2_Status) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) 
dev.off()

summary(lm(SHAP ~ Context, tet2_f, subset = TET2_Status == "Bi-allelic"))
summary(lm(SHAP ~ Context, tet2_f, subset = TET2_Status == "Mono-allelic"))

### ASXL1
asxl1_f <- shap_patient %>%
  filter(Gene == "ASXL1") %>%
  mutate(HR = 2**SHAP,
    Context = ifelse(TP53 == 1, "ASXL1 + TP53",
      ifelse(SF3B1 == 1, "ASXL1 + SF3B1",
        ifelse(RUNX1 == 1, "ASXL1 + RUNX1",
          ifelse(CBL == 1, "ASXL1 + CBL",
          ifelse(U2AF1 == 1, "ASXL1 + U2AF1",
            ifelse(SRSF2 == 1, "ASXL1 + SRSF2",
              ifelse(DNMT3A == 1, "ASXL1 + DNMT3A",
              ifelse(ETNK1 == 1, "ASXL1 + ETNK1",
                ifelse(PHF6 == 1, "ASXL1 + PHF6",
                  ifelse(PTPN11 == 1, "ASXL1 + PTPN11",
                    ifelse(STAG2 == 1, "ASXL1 + STAG2",
                      ifelse(del5q == 1, "ASXL1 + del5q",
                        ifelse(del7 == 1, "ASXL1 + del7",
                          ifelse(delY == 1, "ASXL1 + delY",
                ifelse(N_mutations == 1, "Only ASXL1", 
                  ifelse(N_mutations == 2, "ASXL1 and 1 mut.", "ASXL1 and 2+ mut."))))))))))))))))) %>%
    rowwise() %>%
    mutate(main_mut = isMain("U2AF1", ID, maf))
summary(lm(SHAP ~ ., select(asxl1_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ TP53 + SF3B1 + RUNX1 + CBL + U2AF1 + SRSF2 + DNMT3A + ETNK1 + PHF6 + PTPN11 + STAG2  + del5q + del7 + delY + N_mutations + main_mut, asxl1_f))

png("figures/gnn/SHAP_full_ASXL1.png", width = 2500, height = 1500, res = 300)
ggplot(asxl1_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for ASXL1", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*RUNX1,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*del7,  data = clinical)

### SF3B1
sf3b1_f <- shap_patient %>%
  filter(Gene == "SF3B1") %>%
  mutate(HR = 2**SHAP,
    Context = ifelse(TP53 == 1, "SF3B1 + TP53",
      ifelse(PPM1D == 1, "SF3B1 + PPM1D",
          ifelse(STAG2 == 1, "SF3B1 + STAG2",
                ifelse(del5q == 1, "SF3B1 + del5q",
                  ifelse(N_mutations == 1, "Only SF3B1",
                    ifelse(N_mutations == 2, "SF3B1 and 1 mut.", "SF3B1 and 2+ mut."))))))) %>%
    rowwise() %>%
    mutate(main_mut = isMain("SF3B1", ID, maf))
summary(lm(SHAP ~ ., select(sf3b1_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ TP53 +  PPM1D + STAG2 + del5q, sf3b1_f))

png("figures/gnn/SHAP_full_SF3B1.png", width = 2500, height = 1500, res = 300)
ggplot(sf3b1_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for SF3B1", y = "HR", x = "Context") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ SF3B1*del5q,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ SF3B1*TP53,  data = clinical)

summary(lm(SHAP ~ Context, sf3b1_f))

### DNMT3A
dnmt3a_f <- shap_patient %>%
  filter(Gene == "DNMT3A") %>%
  mutate(HR = 2**SHAP,
    Context = ifelse(TP53 == 1, "DNMT3A + TP53",
      ifelse(SF3B1 == 1, "DNMT3A + SF3B1",
        ifelse(CBL == 1, "DNMT3A + CBL",
          ifelse(U2AF1 == 1, "DNMT3A + U2AF1",
            ifelse(CEBPA == 1, "DNMT3A + CEBPA",
              ifelse(GNB1 == 1, "DNMT3A + GNB1",
                ifelse(NF1 == 1, "DNMT3A + NF1",
          ifelse(N_mutations == 1, "Only DNMT3A",
            ifelse(N_mutations == 2, "DNMT3A and 1 mut.", "DNMT3A and 2+ mut.")))))))))) %>%
    rowwise() %>%
    mutate(main_mut = isMain("DNMT3A", ID, maf))
summary(lm(SHAP ~ ., select(dnmt3a_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ TP53 + SF3B1 + CBL + U2AF1  + CEBPA + GNB1 + NF1 + TET2bi + N_mutations + main_mut, dnmt3a_f))

png("figures/gnn/SHAP_full_DNMT3A.png", width = 2500, height = 1500, res = 300)
ggplot(dnmt3a_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for DNMT3A", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ main_mut, labeller = label_both)
dev.off()

summary(lm(SHAP ~ Context, dnmt3a_f))

coxph(Surv(LFS_YEARS, LFS_STATUS) ~ DNMT3A*U2AF1,  data = clinical)

### SRSF2
srsf2_f <- shap_patient %>%
  filter(Gene == "SRSF2") %>%
  mutate(HR = 2**SHAP,
    Context = case_when(
      RUNX1 == 1 ~ "SRSF2 + RUNX1",
      ASXL1 == 1 ~ "SRSF2 + ASXL1",
      IDH1 == 1 ~ "SRSF2 + IDH1",
      TET2bi == 1 ~ "SRSF2 + TET2 bi",
      TET2other == 1 ~ "SRSF2 + TET2 mono",
      N_mutations == 1 ~ "Only SRSF2",
      N_mutations == 2 ~ "SRSF2 and 1 mut.",
      TRUE ~ "SRSF2 and 2+ mut."
    )) %>%
    rowwise() %>%
    mutate(main_mut = ifelse(isMain("SRSF2", ID, maf), "Main mut.", "Secondary mut."))
summary(lm(SHAP ~ ., select(srsf2_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ TP53 + RUNX1 + ASXL1 + CEBPA + GNB1 + IDH1 + STAG2 + TET2bi + TET2other  + N_mutations + main_mut, srsf2_f))

png("figures/gnn/SHAP_full_SRSF2.png", width = 2500, height = 1500, res = 300)
ggplot(srsf2_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for SRSF2", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ main_mut)
dev.off()

summary(lm(SHAP ~ Context, srsf2_f))
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ SRSF2*ASXL1,  data = clinical)

### RUNX1
runx1_f <- shap_patient %>%
  filter(Gene == "RUNX1") %>%
  mutate(HR = 2**SHAP,
    Context = case_when(
      EZH2 == 1 ~ "RUNX1 + EZH2",
      U2AF1 == 1 ~ "RUNX1 + U2AF1",
      IDH1 == 1 ~ "RUNX1 + IDH1",
      SETBP1 == 1 ~ "RUNX1 + SETBP1",
      N_mutations == 1 ~ "Only RUNX1",
      N_mutations == 2 ~ "RUNX1 and 1 mut.",
      TRUE ~ "RUNX1 and 2+ mut."
    )) %>%
    rowwise() %>%
    mutate(main_mut = isMain("RUNX1", ID, maf))
summary(lm(SHAP ~ ., select(runx1_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ EZH2 + U2AF1 + SRSF2  + KRAS + IDH1 + SETBP1  + CYTO_IPSSR + N_mutations + main_mut, runx1_f))

png("figures/gnn/SHAP_full_RUNX1.png", width = 2500, height = 1500, res = 300)
ggplot(runx1_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for RUNX1", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

summary(lm(SHAP ~ Context, runx1_f))
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ RUNX1*U2AF1,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ RUNX1*EZH2,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ RUNX1*CYTO_IPSSR,  data = clinical)

png("figures/gnn/SHAP_full_RUNX1_CYTO.png", width = 2500, height = 1500, res = 300)
ggplot(runx1_f, aes(x = CYTO_IPSSR, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for RUNX1", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

### TP53
tp53_f <- shap_patient %>%
  filter(Gene == "TP53") %>%
  mutate(HR = 2**SHAP,
    TP53_Status = ifelse(TP53multi == 1, "Multi-allelic", ifelse(TP53mono == 1, "Mono-allelic", "Other")),
    Context = case_when(
      complex == 1 ~ "TP53 + complex",
      del7 == 1 | del7q == 1 ~ "TP53 + del7/del7q",
      N_mutations == 1 ~ "Only TP53",
      N_mutations == 2 ~ "TP53 and 1 mut.",
      TRUE ~ "TP53 and 2+ mut."
    )) %>%
    rowwise() %>%
    mutate(main_mut = isMain("TP53", ID, maf))
summary(lm(SHAP ~ ., select(tp53_f, -c(ID, TP53multi, TP53mono, Gene, Context, HR, TP53_Status))))
summary(lm(SHAP ~  complex + del7 + del7q + CYTO_IPSSR + N_mutations + main_mut, tp53_f))

png("figures/gnn/SHAP_full_TP53.png", width = 2500, height = 1500, res = 300)
ggplot(tp53_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for TP53", y = "HR", x = "Context") +
  theme_bw() +
  facet_wrap(~ TP53_Status) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

summary(lm(SHAP ~ Context, tp53_f, subset = TP53_Status == "Multi-allelic"))
summary(lm(SHAP ~ Context, tp53_f, subset = TP53_Status == "Mono-allelic"))

png("figures/gnn/SHAP_full_TP53_CYTO.png", width = 2500, height = 1500, res = 300)
ggplot(tp53_f, aes(x = CYTO_IPSSR, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for TP53", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  facet_wrap(~ TP53_Status) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

### STAG2
stag2_f <- shap_patient %>%
  filter(Gene == "STAG2") %>%
  mutate(HR = 2**SHAP,
    Context = case_when(
      RUNX1 == 1 ~ "STAG2 + RUNX1",
      NRAS == 1 ~ "STAG2 + NRAS",
      U2AF1 == 1 ~ "STAG2 + U2AF1",
      ASXL1 == 1 ~ "STAG2 + ASXL1",
      BCOR == 1 ~ "STAG2 + BCOR",
      CEBPA == 1 ~ "STAG2 + CEBPA",
      PRPF8 == 1 ~ "STAG2 + PRPF8",
      N_mutations == 1 ~ "Only STAG2",
      N_mutations == 2 ~ "STAG2 and 1 mut.",
      TRUE ~ "STAG2 and 2+ mut."
    )) %>%
    rowwise() %>%
    mutate(main_mut = isMain("STAG2", ID, maf))
summary(lm(SHAP ~ ., select(stag2_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ RUNX1 + NRAS +  U2AF1 + ASXL1 + BCOR + CEBPA + PRPF8 + N_mutations + main_mut, stag2_f))

png("figures/gnn/SHAP_full_STAG2.png", width = 2500, height = 1500, res = 300)
ggplot(stag2_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for STAG2", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

summary(lm(SHAP ~ Context, stag2_f))
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ STAG2*RUNX1,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ STAG2*NRAS,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ STAG2*U2AF1,  data = clinical)


### U2AF1
u2af1_f <- shap_patient %>%
  filter(Gene == "U2AF1") %>%
  mutate(HR = 2**SHAP,
    Context = case_when(
      KRAS == 1 ~ "U2AF1 + KRAS",
      PTPN11 == 1 ~ "U2AF1 + PTPN11",
      del5q == 1 ~ "U2AF1 + del5q",
      N_mutations == 1 ~ "Only U2AF1",
      N_mutations == 2 ~ "U2AF1 and 1 mut.",
      N_mutations == 3 ~ "U2AF1 and 2 mut.",
      TRUE ~ "U2AF1 and 3+ mut."
    )) %>%
    rowwise() %>%
    mutate(main_mut = isMain("U2AF1", ID, maf))
summary(lm(SHAP ~ ., select(u2af1_f, -c(ID, Gene, Context, HR))))
summary(lm(SHAP ~ KRAS  + PTPN11 + del5q  + N_mutations + main_mut, u2af1_f))

png("figures/gnn/SHAP_full_U2AF1.png", width = 2500, height = 1500, res = 300)
ggplot(u2af1_f, aes(x = Context, y = HR)) +
  geom_boxplot() +
  labs(title = "HR for U2AF1", y = "HR", x = "Context") +
  theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

summary(lm(SHAP ~ Context, u2af1_f))
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ U2AF1*del5q,  data = clinical)
| Variable      | TET2 | ASXL1 | SF3B1 | DNMT3A | SRSF2 | RUNX1 | TP53 | STAG2 | U2AF1 |
| ------------- | :--: | :---: | :---: | :----: | :---: | :---: | :--: | :---: | :---: |
| **TET2bi**    |      |       |       |    X   |   X   |       |      |       |       |
| **TET2other** |   X  |       |       |        |   X   |       |      |       |       |
| **ASXL1**     |      |       |       |        |   X   |       |      |   X   |       |
| **SF3B1**     |   X  |   X   |       |    X   |       |       |      |       |       |
| **DNMT3A**    |      |   X   |       |        |       |       |      |       |       |
| **SRSF2**     |   X  |   X   |       |        |       |   X   |      |       |       |
| **RUNX1**     |   X  |   X   |       |        |   X   |       |      |   X   |       |
| **TP53**      |      |   X   |   X   |    X   |   X   |       |      |       |       |
| **STAG2**     |   X  |   X   |   X   |        |   X   |       |      |       |       |
| **U2AF1**     |      |   X   |       |    X   |       |   X   |      |   X   |       |
| CBL           |   X  |   X   |       |    X   |       |       |      |       |       |
| complex       |   X  |       |       |        |       |       |   X  |       |       |
| del7          |   X  |   X   |       |        |       |       |   X  |       |       |
| del7q         |   X  |       |       |        |       |       |   X  |       |       |
| N\_mutations  |   X  |   X   |       |    X   |   X   |   X   |   X  |   X   |   X   |
| main\_mut     |   X  |   X   |       |    X   |   X   |   X   |   X  |   X   |   X   |
| ETNK1         |      |   X   |       |        |       |       |      |       |       |
| PHF6          |      |   X   |       |        |       |       |      |       |       |
| PTPN11        |      |   X   |       |        |       |       |      |       |   X   |
| del5q         |      |   X   |   X   |        |       |       |      |       |   X   |
| delY          |      |   X   |       |        |       |       |      |       |       |
| PPM1D         |      |       |   X   |        |       |       |      |       |       |
| CEBPA         |      |       |       |    X   |   X   |       |      |   X   |       |
| GNB1          |      |       |       |    X   |   X   |       |      |       |       |
| NF1           |      |       |       |    X   |       |       |      |       |       |
| IDH1          |      |       |       |        |   X   |   X   |      |       |       |
| EZH2          |      |       |       |        |       |   X   |      |       |       |
| KRAS          |      |       |       |        |       |   X   |      |       |   X   |
| SETBP1        |      |       |       |        |       |   X   |      |       |       |
| CYTO\_IPSSR   |      |       |       |        |       |   X   |   X  |       |       |
| NRAS          |      |       |       |        |       |       |      |   X   |       |
| BCOR          |      |       |       |        |       |       |      |   X   |       |
| PRPF8         |      |       |       |        |       |       |      |   X   |       |

shap_combs_sel <- shap_patient %>%
  filter(Gene %in% c("TET2","ASXL1","SF3B1","DNMT3A","SRSF2","RUNX1","TP53","STAG2","U2AF1")) %>%
  mutate(
    DNMT3A_TET2bi    = as.integer(DNMT3A    == 1 & TET2bi == 1),
    SRSF2_TET2bi     = as.integer(SRSF2     == 1 & TET2bi == 1),
    SRSF2_TET2other  = as.integer(SRSF2     == 1 & TET2other == 1),
    ASXL1_SRSF2      = as.integer(ASXL1     == 1 & SRSF2 == 1),
    ASXL1_STAG2      = as.integer(ASXL1     == 1 & STAG2 == 1),
    SF3B1_TET2bi     = as.integer(SF3B1     == 1 & TET2bi == 1),
    SF3B1_ASXL1      = as.integer(SF3B1     == 1 & ASXL1 == 1),
    SF3B1_DNMT3A     = as.integer(SF3B1     == 1 & DNMT3A == 1),
    DNMT3A_ASXL1     = as.integer(DNMT3A    == 1 & ASXL1 == 1),
    SRSF2_ASXL1      = as.integer(SRSF2     == 1 & ASXL1 == 1),
    SRSF2_RUNX1      = as.integer(SRSF2     == 1 & RUNX1 == 1),
    RUNX1_TET2bi     = as.integer(RUNX1     == 1 & TET2bi == 1),
    RUNX1_ASXL1      = as.integer(RUNX1     == 1 & ASXL1 == 1),
    RUNX1_SRSF2      = as.integer(RUNX1     == 1 & SRSF2 == 1),
    RUNX1_STAG2      = as.integer(RUNX1     == 1 & STAG2 == 1),
    TP53_ASXL1       = as.integer(TP53      == 1 & ASXL1 == 1),
    TP53_SF3B1       = as.integer(TP53      == 1 & SF3B1 == 1),
    TP53_DNMT3A      = as.integer(TP53      == 1 & DNMT3A == 1),
    TP53_SRSF2       = as.integer(TP53      == 1 & SRSF2 == 1),
    STAG2_TET2bi     = as.integer(STAG2     == 1 & TET2bi == 1),
    STAG2_ASXL1      = as.integer(STAG2     == 1 & ASXL1 == 1),
    STAG2_SF3B1      = as.integer(STAG2     == 1 & SF3B1 == 1),
    STAG2_SRSF2      = as.integer(STAG2     == 1 & SRSF2 == 1),
    U2AF1_ASXL1      = as.integer(U2AF1     == 1 & ASXL1 == 1),
    U2AF1_DNMT3A     = as.integer(U2AF1     == 1 & DNMT3A == 1),
    U2AF1_RUNX1      = as.integer(U2AF1     == 1 & RUNX1 == 1),
    U2AF1_STAG2      = as.integer(U2AF1     == 1 & STAG2 == 1)
  )


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

save(shaps, shaps_per_gene, comb_shap, file = "results/gnn/v1_boolean_gene_SHAP.Rdata")

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

png("figures/gnn/gene_embedding_heatmap.png", width = 700, height = 1200)
pheatmap::pheatmap(filt_mat)
dev.off()

filt_mat2 <- filt_mat[unique(comb_shap$Gene), ]  
png("figures/gnn/gene_embedding_heatmap_gene_mod.png", width = 700, height = 1000)
pheatmap::pheatmap(filt_mat2)
dev.off()

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


## Explore further failing genes
load("results/preprocess/clinical_preproc.Rdata")

comb_shap_patient <- left_join(comb_shap,
  clinical %>% select(ID, any_of(IPSSM_genes), TP53multi, TP53mono, 
    TET2bi, TET2other, N_mutations, complex, del5q, plus8, del20q, del7, del7q, delY,
    CYTO_IPSSR), by = "ID")

## Explore TP53
tp53 <- filter(comb_shap_patient, Gene == "TP53")   %>%
  mutate(TP53_Status = ifelse(TP53multi == 1, "Multi", ifelse(TP53mono == 1, "Mono", "Other")),
        Mutation_status = 
        ifelse(del7 == 1, "TP53 + del7",
          ifelse(del7q == 1, "TP53 + del7q",
          ifelse(ASXL1 == 1, "TP53 + ASXL1",
                    ifelse(N_mutations == 1, "Only TP53",
                    ifelse(N_mutations == 2, "TP53 and 1 mut.",  "TP53 and 2+ mut."))))),
        Mutation_status = factor(Mutation_status, levels = c("Only TP53",  "TP53 + del7", 
        "TP53 + del7q", "TP53 + ASXL1", "TP53 and 1 mut.", 
        "TP53 and 2+ mut."))
  )

summary(lm(SHAP_diff ~ ., select(tp53, -c(ID, TP53mono, TP53, Gene, SHAP_gene, Mutation_status, TP53_Status, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~  DNMT3A + ASXL1 +   complex +  del7 + del7q + CYTO_IPSSR, tp53))

png("figures/gnn/SHAP_TP53.png", width = 2000, height = 2000, res = 300)
tp53 %>%
 ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for TP53", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  facet_grid( ~ TP53_Status, scales = "free") +
  theme_bw() 
dev.off()

tp53 %>%
  group_by(TP53_Status, Mutation_status) %>%
  summarize(cor = cor(SHAP_gene, SHAP_full, use = "complete.obs")) 

png("figures/gnn/SHAP_TP53_complex.png", width = 2000, height = 2000, res = 300)
tp53 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = complex)) +
  geom_point() +
  labs(title = "SHAP Value Differences for TP53", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  facet_grid( ~ TP53_Status, scales = "free") +
  theme_bw() 
dev.off()

coxph(Surv(LFS_YEARS, LFS_STATUS) ~ TP53multi*ASXL1,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ TP53multi*N_mutations,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ TP53mono*ASXL1,  data = clinical)

png("figures/gnn/SHAP_TP53_CYTO.png", width = 2000, height = 2000, res = 300)
tp53 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = CYTO_IPSSR)) +
  geom_point() +
  labs(title = "SHAP Value Differences for TP53", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  facet_grid( ~ TP53_Status, scales = "free") +
  theme_bw() 
dev.off()



## Explore SF3B1
sf3b1 <- filter(comb_shap_patient, Gene == "SF3B1")  %>%
  mutate(Mutation_status = 
    ifelse(del5q == 1, "SF3B1 + del5q",
          ifelse(STAG2 == 1, "SF3B1 + STAG2",
            ifelse(RUNX1 == 1, "SF3B1 + RUNX1", 
              ifelse(SRSF2 == 1, "SF3B1 + SRSF2",
              ifelse(TET2bi == 1, "SF3B1 + TET2bi",
                ifelse(IDH2 == 1, "SF3B1 + IDH2",
                ifelse(EZH2 == 1, "SF3B1 + EZH2",
                ifelse(TET2other == 1, "SF3B1 + TET2other",
                ifelse(N_mutations == 1, "Only SF3B1", 
                  ifelse(N_mutations == 2, "SF3B1 and 1 mut.", "SF3B1 and 2+ mut.")))))))))),
        Mutation_status = factor(Mutation_status, 
          levels = c("Only SF3B1", "SF3B1 + del5q", "SF3B1 + STAG2", "SF3B1 + RUNX1", 
               "SF3B1 + TET2bi", "SF3B1 + IDH2", "SF3B1 + SRSF2", "SF3B1 + EZH2", 
               "SF3B1 + TET2other", "SF3B1 and 1 mut.", "SF3B1 and 2+ mut.")))
summary(lm(SHAP_diff ~ ., select(sf3b1, -c(ID, TP53mono, SF3B1, Gene, SHAP_gene, Mutation_status, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ RUNX1 + IDH2 + EZH2 + SRSF2  + STAG2 + TET2bi + TET2other  + del5q  + CYTO_IPSSR, sf3b1))

png("figures/gnn/SHAP_SF3B1.png", width = 3000, height = 2200, res = 300)
sf3b1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for SF3B1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  facet_wrap(~ Mutation_status, scales = "free") 
dev.off()


## Explore DNMT3A
dnmt3a <- filter(comb_shap_patient, Gene == "DNMT3A")  %>%
  mutate(Mutation_status = 
  ifelse(EZH2 == 1, "DNMT3A + EZH2",
    ifelse(STAG2 == 1, "DNMT3A + STAG2",
      ifelse(TET2other == 1, "DNMT3A + TET2other",
            ifelse(N_mutations == 1, "Only DNMT3A", 
              ifelse(N_mutations == 2, "DNMT3A and 1 mut.", "DNMT3A and 2+ mut."))))),
        Mutation_status = factor(Mutation_status, 
          levels = c("Only DNMT3A", "DNMT3A + EZH2", "DNMT3A + STAG2", "DNMT3A + TET2other", "DNMT3A and 1 mut.", "DNMT3A and 2+ mut.")))
summary(lm(SHAP_diff ~ ., select(dnmt3a, -c(ID, TP53mono,Mutation_status, DNMT3A, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ EZH2 + STAG2 + TET2other + N_mutations + CYTO_IPSSR, dnmt3a))

png("figures/gnn/SHAP_DNMT3A.png", width = 3000, height = 2200, res = 300)
dnmt3a %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for DNMT3A", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw()
dev.off()

png("figures/gnn/SHAP_DNMT3A_CYTO.png", width = 3000, height = 2200, res = 300)
dnmt3a %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = CYTO_IPSSR)) +
  geom_point() +
  labs(title = "SHAP Value Differences for DNMT3A", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw()
dev.off()

## Explore TET2
tet2 <- filter(comb_shap_patient, Gene == "TET2")  %>%
  mutate(TET2_Status = ifelse(TET2bi == 1, "TET2bi", ifelse(TET2other == 1, "TET2other", "Other")),
         Mutation_status = 
           ifelse(NF1 == 1, "TET2 + NF1",
            ifelse(RUNX1 == 1, "TET2 + RUNX1",
                    ifelse(N_mutations == 1, "Only TET2", 
                      ifelse(N_mutations == 2, "TET2 and 1 mut.", "TET2 and 2+ mut.")))),
         Mutation_status = factor(Mutation_status,
           levels = c("Only TET2", "TET2 + NF1", "TET2 + RUNX1", "TET2 and 1 mut.", "TET2 and 2+ mut.")))
summary(lm(SHAP_diff ~ ., select(tet2, -c(ID, TP53mono, TET2bi, Mutation_status, TET2_Status, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ RUNX1 + NF1  + N_mutations +  del7q, tet2))

png("figures/gnn/SHAP_TET2.png", width = 3000, height = 2200, res = 300)
tet2 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for TET2", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  facet_grid( ~ TET2_Status)
dev.off()

## U2AF1


u2af1 <- filter(comb_shap_patient, Gene == "U2AF1")  %>%
  mutate(Mutation_status = 
    ifelse(ETNK1 == 1, "U2AF1 + ETNK1",
      ifelse(STAG2 == 1, "U2AF1 + STAG2",
        ifelse(complex == "complex", "U2AF1 + complex",
          ifelse(N_mutations == 1, "Only U2AF1", 
            ifelse(N_mutations == 2, "U2AF1 and 1 mut.", 
              ifelse(N_mutations == 3, "U2AF1 and 2 mut.", "U2AF1 and 3+ mut.")))))),
    Mutation_status = factor(Mutation_status,
      levels = c("Only U2AF1", "U2AF1 + ETNK1", "U2AF1 + STAG2", "U2AF1 + complex", "U2AF1 and 1 mut.", "U2AF1 and 2 mut.", "U2AF1 and 3+ mut."))) %>%
    rowwise() %>%
    mutate(main_mut = isMain("U2AF1", ID, maf)) 
summary(lm(SHAP_diff ~ ., select(u2af1, -c(ID, TP53mono, Mutation_status, U2AF1, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~   ETNK1 + STAG2   + complex + N_mutations, u2af1))

png("figures/gnn/SHAP_U2AF1.png", width = 3000, height = 2200, res = 300)
u2af1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for U2AF1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  facet_wrap(~ Mutation_status)
dev.off()

png("figures/gnn/SHAP_U2AF1b.png", width = 3000, height = 2200, res = 300)
u2af1 %>%
  mutate(Status2 = ifelse(Mutation_status %in% c("U2AF1 + ETNK1", "U2AF1 + complex", "U2AF1 + STAG2"), "Def Comb.", as.character(Mutation_status))) %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Status2)) +
  geom_point() +
  labs(title = "SHAP Value Differences for U2AF1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() 
dev.off()

png("figures/gnn/SHAP_U2AF1_mainmut.png", width = 3000, height = 2200, res = 300)
u2af1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for U2AF1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  facet_grid( ~ main_mut)
dev.off()
u2af1$class <- modelo_gmm$classification
cot <- multinom(factor(class) ~ .,  select(u2af1, -c(ID, TP53mono, Mutation_status, U2AF1, Gene, SHAP_gene, Testing_Gene, SHAP_full, SHAP_diff)))

## ASXL1
asxl1 <- filter(comb_shap_patient, Gene == "ASXL1")  %>%
  mutate(
    Mutation_status = case_when(
      NRAS == 1 ~ "ASXL1 + NRAS",
      ETV6 == 1 ~ "ASXL1 + ETV6",
      EZH2 == 1 ~ "ASXL1 + EZH2",
      SRSF2 == 1 ~ "ASXL1 + SRSF2",
      STAG2 == 1 ~ "ASXL1 + STAG2",
      SF3B1 == 1 ~ "ASXL1 + SF3B1",
      N_mutations == 1 ~ "Only ASXL1",
      N_mutations == 2 ~ "ASXL1 and 1 mut.",
      TRUE ~ "ASXL1 and 2+ mut."
    ),
    Mutation_status = factor(
      Mutation_status,
      levels = c(
        "Only ASXL1",
        "ASXL1 + SF3B1",
        "ASXL1 + NRAS",
        "ASXL1 + ETV6",
        "ASXL1 + EZH2",
        "ASXL1 + SRSF2",
        "ASXL1 + STAG2",
        "ASXL1 and 1 mut.",
        "ASXL1 and 2+ mut."
      )
    )
  )
summary(lm(SHAP_diff ~ ., select(asxl1, -c(ID, TP53mono, Mutation_status, ASXL1, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ SF3B1 + NRAS + ETV6 + EZH2 + SRSF2 + BCOR + ETNK1 + NF1 + STAG2 + N_mutations, asxl1))


png("figures/gnn/SHAP_ASXL1.png", width = 3000, height = 2200, res = 300)
asxl1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for ASXL1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  facet_wrap(~ Mutation_status)
dev.off()
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*SF3B1,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*SRSF2,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*STAG2,  data = clinical)
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ASXL1*EZH2,  data = clinical)

## NF1
nf1 <- filter(comb_shap_patient, Gene == "NF1")  %>%
  mutate(Mutation_status = 
    ifelse(N_mutations == 1, "Only NF1", 
      ifelse(N_mutations == 2, "NF1 and 1 mut.", "NF1 and 2+ mut.")),
    Mutation_status = factor(Mutation_status,
      levels = c("Only NF1", "NF1 and 1 mut.", "NF1 and 2+ mut.")))
summary(lm(SHAP_diff ~ ., select(nf1, -c(ID, TP53mono, Mutation_status, NF1, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ CEBPA + N_mutations, nf1))

png("figures/gnn/SHAP_NF1.png", width = 2000, height = 1600, res = 300)
nf1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for NF1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() 
dev.off()

## RUNX1
runx1 <- filter(comb_shap_patient, Gene == "RUNX1")  %>%
  mutate(Mutation_status = 
  ifelse(SRSF2 == 1, "RUNX1 + SRSF2",
    ifelse(STAG2 == 1, "RUNX1 + STAG2",
    ifelse(N_mutations == 1, "Only RUNX1", 
      ifelse(N_mutations == 2, "RUNX1 and 1 mut.", "RUNX1 and 2+ mut.")))),
    Mutation_status = factor(Mutation_status,
      levels = c("Only RUNX1", "RUNX1 + SRSF2", "RUNX1 + STAG2", "RUNX1 and 1 mut.", "RUNX1 and 2+ mut.")))
summary(lm(SHAP_diff ~ ., select(runx1, -c(ID, TP53mono, Mutation_status, RUNX1, Gene, SHAP_gene, Testing_Gene, SHAP_full))))
summary(lm(SHAP_diff ~ SRSF2  + STAG2 +  CYTO_IPSSR + N_mutations, runx1))

png("figures/gnn/SHAP_RUNX1.png", width = 2000, height = 1600, res = 300)
runx1 %>%
  ggplot(aes(x = SHAP_gene, y = SHAP_full, color = Mutation_status)) +
  geom_point() +
  labs(title = "SHAP Value Differences for RUNX1", y = "SHAP Full", x = "SHAP Gene") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw()
dev.off()

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