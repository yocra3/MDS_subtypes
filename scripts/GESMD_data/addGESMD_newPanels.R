## Check new data
library(readxl)
library(tidyverse)

## Load GESMD data ####
load("results/gesmd_data_all.Rdata")

## Load new panels ####
new_pat <- read_xlsx("data/GESMD_genetics/RESMD Molecular_TODOS_12-11-24.xlsx") %>%
  mutate(register_number  = as.integer(`ID paciente RESMD`),
         vaf = as.double(`VAF (%)`),
         vaf = ifelse(vaf < 1, vaf * 100, vaf), 
         panel_name = `Nombre del panel estudiado (indicar detalles del panel en pestaña 2)`)
panels_info <- read_xlsx("data/GESMD_genetics/RESMD Molecular_TODOS_12-11-24.xlsx", 
                         sheet = "Paneles", .name_repair = "minimal")

panels <- list()
cols <- colnames(panels_info)
for (idx in seq_len(ncol(panels_info))){
  
  col <- cols[idx]
  if (col == "Nombre del panel"){
    panel_name <- unlist(panels_info[1, idx])
  }
  if (col == "Lista de genes"){
    genes <- panels_info[[idx]]
    genes <- genes[!is.na(genes)]
    panels[[panel_name]] <- genes
    
  }
}

## Select all tested genes
all_genes <- unlist(panels) %>% unique()

# Process panels ####

## Select panels
paciente_panel <- select(new_pat, register_number, panel_name) %>%
  distinct() %>%
  mutate(panel_name = case_match(panel_name,
                                 "Panel Oncomine" ~ "Oncomine Myeloid Research_Thermo Fisher",
                                 "Panel Ampliseq" ~ "Ion AmpliSeq AML Research Panel",
                                 .default = panel_name))

## Create empty matrix
mut_mat <- matrix(NA, nrow = nrow(paciente_panel), ncol = length(all_genes),
                  dimnames = list(paciente_panel$register_number, all_genes))

## Fill present genes in panels ####
for (panel in unique(paciente_panel$panel_name)){
  
  panel_patients <- subset(paciente_panel, panel_name == panel)$register_number
  
  sel_genes <- intersect(panels[[panel]], colnames(mut_mat))
  
  mut_mat[as.character(panel_patients), sel_genes] <- 0
}


pat_tab <- select(new_pat, register_number, Gen, vaf) %>%
  group_by(register_number, Gen) %>%
  summarize(mut = sum(vaf)) %>%
  mutate(mut = ifelse(!Gen %in% c("TP53", "TET2") & mut > 1, 1, mut),
         Gen = ifelse(Gen == "TP53+HF165:H196", "TP53", Gen)) %>%
  filter(Gen %in% all_genes & !is.na(Gen))

for (patient in unique(pat_tab$register_number)){
  
  inputs <- subset(pat_tab, register_number == patient)
  
  mut_mat[as.character(patient), inputs$Gen] <- inputs$mut
} 

## Convert to tibble and modify variables
gesmd_mut <- mut_mat %>%
  as_tibble() %>%
  mutate(register_number = rownames(mut_mat),
         TP53mut = TP53 > 0, 
         TP53loh = NA,
         MLL_PTD = NA,
         TET2bi = ifelse(TET2 >= 50, 1, 0),
         TET2other = ifelse(TET2 > 0 & TET2 < 50, 1, 0),
         TP53multi = ifelse(TP53 >= 50, 1, 0),
         TP53mono = ifelse(TP53 > 0 & TP53 < 50, 1, 0),
         TET2 = ifelse(TET2 > 0, 1, 0),
         TP53 = ifelse(TP53 > 0, 1, 0)) 


tp53_tab <- select(new_pat, register_number, Gen, vaf) %>%
  filter(Gen %in% c("TP53", "TP53+HF165:H196")) %>%
  group_by(register_number) %>%
  summarize(n_rel = sum(vaf >= 10),
            vaf = sum(vaf)) %>%
  mutate(TP53multi = ifelse(n_rel >= 2, 1, 0),
         register_number = as.character(register_number)) %>%
  select(register_number, TP53multi)

gesmd_mut <- gesmd_mut %>%
  rows_update(y = tp53_tab, by = "register_number") %>%
  mutate(TP53mono = ifelse(TP53mono == 1 & TP53multi == 1, 0, TP53mono))

save(gesmd_mut, file = "results/gesmd_update_0626.Rdata")

gesmd_data_new <- gesmd_data %>%
  rows_update(
  y = gesmd_mut,
  by = "register_number"
)  %>%
  mutate(consensus = ifelse(is.na(TP53) | is.na(SF3B1), NA,  
           ifelse(TP53multi == 1 & BM_BLAST <= 20, "Mutated TP53",
                            ifelse(del5q == 1 & del7q == 0 & BM_BLAST <= 5, "del5q",
                                   ifelse(SF3B1 > 0 & del7q == 0 & complex == 0 & BM_BLAST <= 5, "mutated SF3B1",
                                          ifelse(BM_BLAST <= 5, "Low blasts",
                                                 ifelse(BM_BLAST > 10, "MDS-IB2",
                                                        ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other"))))))))

