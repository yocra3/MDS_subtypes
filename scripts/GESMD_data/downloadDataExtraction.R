library(tidyverse)
library(readxl)
library(RPostgres)

## Get ids
patients <- read_xlsx("data/datosCUN/Mielodisplasias.xlsx")
ids <- paste0(patients$NH, collapse = ",")

#primero conectar tunel SSH via fichero bash
ssh -i private.pem -N -L 5556:prehm-pro-rds-i.cfeequry4htk.eu-west-3.rds.amazonaws.com:5432 ec2-user@35.181.32.109

## Connect to db
drv <- RPostgreSQL::PostgreSQL()
con <- dbConnect(RPostgres::Postgres(), user="cruizarenas_resmdro",
                 password='s"4/Z7bG=m%9j3W',
                 dbname="masihdas_resmd",
                 port=5556, # el puerto es el mismo al que rediriges el tunel ssh
                 host="localhost")

all_patients <- dbGetQuery(con, sprintf(
  "select * 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  join current_status cs on p.id = cs.id
  join diagnosis_data dd on p.id = dd.id 
  join family_history_data fhd on p.id = fhd.id
  join personal_history_data phd on p.id = phd.id
  join prognostic_scoring ps on p.id = ps.diagnosis_data 
  join sorror_index s on p.id = s.id
  where p.nip = ANY(ARRAY[%s]::varchar[])", ids)) 

pat_ids <- paste0(all_patients$id, collapse = ",")


all_tables_l1 <- dbGetQuery(con, sprintf("
SELECT
table_schema,
table_name,
column_name
FROM information_schema.columns
WHERE column_name = 'patient_id'
ORDER BY table_name;
"))

join_l1 <- lapply(all_tables_l1$table_name, function(tab_name){
  dbGetQuery(con, sprintf(
    "select * 
  from %s tab
  where tab.patient_id = ANY(ARRAY[%s]::bigint[])", tab_name, pat_ids)) 
})
names(join_l1) <- all_tables_l1$table_name
sapply(join_l1, nrow)

analytical_id <- paste0(join_l1$analytical_data$id, collapse = ",")

bio_data <-   dbGetQuery(con, sprintf(
  "select * 
  from biochemistry_data b
  where b.bio_data_id = ANY(ARRAY[%s]::bigint[])", analytical_id))

cyto_data <-   dbGetQuery(con, sprintf(
  "select * 
  from citogenetical_data c
  where c.cito_data_id = ANY(ARRAY[%s]::bigint[])", analytical_id))


hemogram <- dbGetQuery(con, sprintf(
  "select * 
  from hemogram_data hd
  where hd.hemogram_data_id = ANY(ARRAY[%s]::bigint[])", analytical_id))


molecular <- dbGetQuery(con, sprintf(
  "select * 
  from molecular_data md
  join mutation m on  md.molec_data_id = m.molecular_data_id 
  where md.molec_data_id = ANY(ARRAY[%s]::bigint[])", analytical_id))

morphological <- dbGetQuery(con, sprintf(
  "select * 
  from morphological_data mod
  where mod.morphologycal_data_id = ANY(ARRAY[%s]::bigint[])", analytical_id))

treatment_id <- paste0(join_l1$treatment$id, collapse = ",")

treatment <-  dbGetQuery(con, sprintf(
  "select * 
  from treatment_data t
  join treatment_drug_data tcd on t.id = tcd.treatment_data_id
  where t.id = ANY(ARRAY[%s]::bigint[])", treatment_id))

