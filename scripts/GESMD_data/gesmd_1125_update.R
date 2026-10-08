## Update GESMD data ####
library(tidyverse)

load("results/gesmd_data_all.Rdata")
load( "results/gesmd_cun_update_1125.Rdata")
load( "results/gesmd_update_0626.Rdata")

## Update Mutations
com_cols_mut <- intersect(colnames(gesmd_mut), colnames(gesmd_data))

gesmd_data_mut <- gesmd_data %>%
  rows_patch(
    y = select(gesmd_mut, all_of(com_cols_mut)),
    by = "register_number"
  )  

## Update CUN
com_cols_cun <- intersect(colnames(gesmd_cun), colnames(gesmd_data))

gesmd_cun2 <- mutate(gesmd_cun, TP53mono = ifelse(TP53 == 1 & TP53multi == 0, 1, 0))

gesmd_data_0626 <- gesmd_data_mut %>%
  rows_update(
    y = select(gesmd_cun2, all_of(com_cols_cun)),
    by = "register_number"
  )  %>%
  mutate(consensus = ifelse(is.na(TP53) | is.na(SF3B1) | is.na(del5q), NA,  
                            ifelse((TP53multi == 1 | (TP53mono == 1 & del17p13 == 1)) & BM_BLAST <= 20, "MDS-TP53",
                                   ifelse(del5q == 1 & del7q == 0 & del7 == 0 & BM_BLAST <= 5 & N_aberrations < 3 & complex == 0, "MDS-del5q",
                                          ifelse(SF3B1 == 1 & del7q == 0 & del7 == 0 & RUNX1 == 0 & complex == 0 & BM_BLAST <= 5, "MDS-SF3B1",
                                                 ifelse(BM_BLAST <= 5, "MDS-LB",
                                                        ifelse(BM_BLAST > 10, "MDS-IB2",
                                                               ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other"))))))))
                            

## Eliminar muestras del ICO en el IWS
ICO_IWS <- read_xlsx("data/GESMD_genetics/Muetsras HUGTIP enviadas IPSS-M.xlsx")
ICO_NHC <- ICO_IWS$NHC

drv <- RPostgreSQL::PostgreSQL()
con <- dbConnect(RPostgres::Postgres(), user="cruizarenas_resmdro",
                 password='s"4/Z7bG=m%9j3W',
                 dbname="masihdas_resmd",
                 port=5556, # el puerto es el mismo al que rediriges el tunel ssh
                 host="localhost")
all_patients_match <- dbGetQuery(con,
  "select register_number, nip, h.name as hospital 
  from patient p 
  join hospital h on p.hospital_id = h.id") %>%
  as_tibble()

rep_register <- subset(all_patients_match, nip %in% ICO_NHC & hospital == "ICO Badalona")

gesmd_data_0626_filt <- subset(gesmd_data_0626,!register_number %in% rep_register$register_number)
save(gesmd_data_0626_filt, file = "results/gesmd_data_0626.Rdata")
