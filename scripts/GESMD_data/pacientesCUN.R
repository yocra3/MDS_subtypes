library(tidyverse)
library(lubridate)
library(readxl)
library(RPostgres)
library(xlsx)

#primero conectar tunel SSH via fichero bash
drv <- RPostgreSQL::PostgreSQL()
con <- dbConnect(RPostgres::Postgres(), user="cruizarenas_resmdro",
                 password='s"4/Z7bG=m%9j3W',
                 dbname="masihdas_resmd",
                 port=5556, # el puerto es el mismo al que rediriges el tunel ssh
                 host="localhost")

all_patients <- dbGetQuery(con, "SELECT id, register_number FROM patient")


all_register <- paste0(all_patients$register_number, collapse = ",")

all_person <- dbGetQuery(con, sprintf(
  "select register_number, nip, h.name as hospital, li.code as gender, birthdate 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  where register_number = any(array[%s])", all_register)) %>%
  mutate(SEX = ifelse(gender == "SEXO_HOMBRE", "M", "F")) %>%
  as_tibble()

all_cun <- subset(all_person , hospital == "Clinica Universidad de Navarra") 

cima_labs <- read_xlsx("data/datosCUN/Panel NGS casos SMD CUN 20260616.xlsx") %>%
  mutate(nip = `Nº Historia`,
         nip = gsub("-CUN", "", nip))


cima_labs_com <- subset(cima_labs, nip %in% all_cun$nip)

write.table(all_cun, file = "CUN_patients_GESMD.txt", row.names = FALSE, quote = FALSE)
write.table(cima_labs_com, file = "CIMA_LABS_CUN_patients_GESMD.txt", row.names = FALSE, quote = FALSE)
write.table(cima_labs_com, file = "CIMA_LABS_CUN_patients_GESMD.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)

write.xlsx(as.data.frame(cima_labs_com[, colnames(cima_labs_com) != "Paciente"]),
                 file = "CIMA_LABS_CUN_patients_GESMD.xlsx", row.names = FALSE, showNA = FALSE)
