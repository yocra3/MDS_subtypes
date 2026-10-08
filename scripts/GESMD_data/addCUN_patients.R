## Add CUN missing data ####
library(readxl)
library(tidyverse)

## Load GESMD data
load("results/gesmd_data_all.Rdata")

## Load NIC ####
all_person <- dbGetQuery(con, sprintf(
  "select register_number, nip, h.name as hospital, li.code as gender, birthdate 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  where register_number = any(array[%s])", all_register)) %>%
  mutate(SEX = ifelse(gender == "SEXO_HOMBRE", "M", "F")) %>%
  as_tibble()

cun <- subset(all_person, hospital == "Clinica Universidad de Navarra") %>%
  mutate(nip = gsub(" ", "", nip))


## CUN data ####
## Map to register
cun_pat1 <- read_xlsx("data/GESMD_genetics/Datos SMD 19112025.xlsx") %>%
  mutate(nip = `Número de Historia`,
         nip = gsub("-CUN", "", nip),
         nip = gsub(" ", "", nip),
         momento = `Datos muestra`,
         momento = ifelse(is.na(momento), "Al diagnóstico", momento))

cun_pat_filt <- filter(cun_pat1, Muestra %in% c("DNA", "Médula ósea") & !is.na(nip)) %>%
  filter(momento %in% c("Al diagnóstico", "Al diagnóstico bis")) %>%
  group_by(nip) %>%
  arrange(momento) %>%
  slice_head(n = 1) %>%
  ungroup()

cun_pat_register <- inner_join(cun_pat_filt, cun, by = "nip") 

## UMBRELLA ####
cun_pat2 <- read_xlsx("data/GESMD_genetics/Datos SMD-UMBRELLA 19112025.xlsx",
                      col_types = c(rep("guess", 13), "date",  rep("guess", 3))) %>%
  mutate(nip = `Número de Historia`,
         nip = gsub("-CUN", "", nip),
         nip = gsub(" ", "", nip),
         momento = `Datos muestra`,
         momento = ifelse(is.na(momento), "Al diagnóstico", momento), 
         UMBRELLA_ID = sapply(strsplit(`Pacientes::Nombre`, ","), `[`, 1),
         UMBRELLA_ID = gsub(" UMBRELLA", "", UMBRELLA_ID))
mapping <-  read_xlsx("data/GESMD_genetics/UMBRELLA RESMD FACTURACION 23102025.xlsx") %>%
  mutate(register_number = gsub("RESMD_", "", `Cód./Ref. Interna`),
         UMBRELLA_ID = Apellidos) %>%
  select(UMBRELLA_ID, register_number )

umbrella_register <- left_join(cun_pat2, mapping, by = "UMBRELLA_ID") %>%
  filter(momento %in% c("Al diagnóstico", "Al diagnóstico bis") & Muestra %in% c("DNA", "Médula ósea")) %>%
  group_by(register_number) %>%
  arrange(momento) %>%
  slice_head(n = 1) %>%
  ungroup()


## Algunos errores probables en la BBDD

## Create gene information ####
panel_genes <- list(
  SOPHIA = c(
    "ANKRD26", "ASXL1", "ATRX", "BCOR", "BCORL1", "CALR", 
    "CBL", "CEBPA", "CSF3R", "CSNK1A1", "CUX1", "DDX41", 
    "DNMT3A", "ETNK1", "ETV6", "EZH2", "FLT3", "GATA1", 
    "GATA2", "IDH1", "IDH2", "IKZF1", "JAK2", "KIT", 
    "KMT2A", "KRAS", "MPL", "NF1", "NPM1", "NRAS", 
    "PHF6", "PPM1D", "PTPN11", "RAD21", "RUNX1", "SETBP1", 
    "SF3B1", "SH2B3/LNK", "SMC1A", "SMC3", "SRP72", 
    "SRSF2", "STAG2", "TET2", "TP53", "U2AF1", "WT1", 
    "ZRSR2"
  ), 
  SOPHIAv2 = c(
    "ACD", "ANKRD26", "ASXL1", "ATRX", "BCOR", "BCORL1", 
    "CALR", "CBL", "CEBPA", "CSF3R", "CSNK1A1", "CUX1", 
    "DDX41", "DHX34", "DNMT3A", "ETNK1", "ETV6", "EZH2", 
    "FLT3", "GATA1", "GATA2", "IDH1", "IDH2", "IKZF1", 
    "JAK2", "KIT", "KMT2A", "KRAS", "MBD4", "MECOM", 
    "MPL", "NF1", "NPM1", "NRAS", "PHF6", "PPM1D", 
    "PTPN11", "RAD21", "RUNX1", "SAMD9", "SAMD9L", 
    "SETBP1", "SF3B1", "SH2B3/LNK", "SMC1A", "SMC3", 
    "SRP72", "SRSF2", "STAG2", "TERC", "TERT", "TET2", 
    "TP53", "U2AF1", "WT1", "ZRSR2"
  ),
  TWIST = c(
    "ABL1", "ANKRD26", "APC", "ASXL1", "ATG2B", "ATOX1", "BCOR", "BCORL1", 
    "BRAF", "CALR", "CAV1", "CBL", "CDKN2A", "CDKN2B", "CEBPA", "CREBBP", 
    "CSF2RA", "CSF3R", "CUX1", "CYP11B2", "CYP3A5", "DDX41", "DNMT3A", "EGR1", 
    "ERCC6L2", "ERG", "ETNK1", "ETV6", "EZH2", "FLT3", "GATA1", "GATA2", 
    "GNB1", "GSKIP", "HNF4A", "IDH1", "IDH2", "IKZF1", "IKZF2", "IKZF3", 
    "IL2RB", "IL3RA", "IL7R", "JAK1", "JAK2", "JAK3", "KIT", "KMT2A", "KRAS", 
    "LIG4", "MPL", "MYC", "NBN", "NF1", "NIPBL", "NOTCH1", "NPM1", "NRAS", 
    "NRG1", "P2RY8", "PAX5", "PHF6", "PIGA", "PPM1D", "PRDM14", "PRPF8", 
    "PTEN", "PTPN11", "RAD21", "RRAS", "RUNX1", "SAMD9", "SAMD9L", "SETBP1", 
    "SF3B1", "SH2B3", "SMC1A", "SMC3", "SRC", "SRP72", "SRSF2", "STAG2", 
    "TCOF1", "TERC", "TERT", "TET2", "TP53", "U2AF1", "UBA1", "WT1", "XPC", 
    "ZEB2", "ZRSR2"
  )

)
all_genes <- unlist(panel_genes) %>% unique()
# sel_muts <- c("TET2", "ASXL1", "SRSF2", "DNMT3A", "RUNX1", 
#               "STAG2", "U2AF1", "EZH2", "ZRSR2",  
#               "SF3B1", "TP53", "FLT3", "BCOR", "BCORL1", "CBL",
#               "CEBPA",	"ETV6",	"IDH1",	"IDH2", "KRAS", "NF1", "NPM1",
#               "NRAS",	"SETBP1",	"ETNK1", "GATA2",	"GNB1",	"PHF6",	"PPM1D", "PRPF8",
#               "PTPN11",	"WT1", "MLL_PTD")

com_columns <- intersect(colnames(cun_pat_register), colnames(umbrella_register))

cun_all <- rbind(cun_pat_register[, com_columns], umbrella_register[, com_columns]) %>%
  mutate(panel = ifelse(grepl("Panel NGS PMPv2-SOPHIA", Observaciones), "SOPHIAv2",
                      ifelse(grepl("Panel NGS PANMieloide-SOPHIA", Observaciones), "SOPHIA", 
                             ifelse(grepl("NGS TWIST", Observaciones), "TWIST", NA))),
         register_number = as.character(register_number)) %>%
  filter(!is.na(register_number)) %>%
  group_by(register_number) %>%
  arrange(panel) %>%
  slice_tail(n = 1) %>%
  ungroup()



mut_mat <- matrix(NA, nrow = nrow(cun_all), ncol = length(all_genes),
                  dimnames = list(as.character(cun_all$register_number), all_genes))

## Fill present genes in panels ####
for (panel in unique(cun_all$panel)){
  
  panel_patients <- subset(cun_all, panel == panel)$register_number
  
  sel_genes <- intersect(panel_genes[[panel]], colnames(mut_mat))
  
  mut_mat[as.character(panel_patients), sel_genes] <- 0
}

gesmd_cun <- as_tibble(mut_mat) %>%
  mutate(register_number = rownames(mut_mat),
  across(
    .cols = all_of(all_genes),
    .fns = ~ifelse(grepl(cur_column(), cun_all$Observaciones), 1, 0)
  ),
  TET2bi = ifelse(grepl("TET2 (2)", cun_all$Observaciones, fixed = TRUE), 1, 0),
  TET2other = ifelse(TET2 == 1 & TET2bi == 0, 1, 0),
  TP53multi = ifelse(grepl("TP53 (2)", cun_all$Observaciones, fixed = TRUE), 1, 0)  )

save(gesmd_cun, file = "results/gesmd_cun_update_1125.Rdata")
