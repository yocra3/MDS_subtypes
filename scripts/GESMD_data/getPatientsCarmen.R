library(tidyverse)
library(lubridate)
library(readxl)
#library(ssh)
library(RPostgres)

#primero conectar tunel SSH via fichero bash
ssh -i private.pem -N -L 5556:prehm-pro-rds-i.cfeequry4htk.eu-west-3.rds.amazonaws.com:5432 ec2-user@35.181.32.109

drv <- RPostgreSQL::PostgreSQL()
con <- dbConnect(RPostgres::Postgres(), user="cruizarenas_resmdro",
                 password='s"4/Z7bG=m%9j3W',
                 dbname="masihdas_resmd",
                 port=5556, # el puerto es el mismo al que rediriges el tunel ssh
                 host="localhost")

## Get patients ####
all_patients <- dbGetQuery(con, "SELECT id, register_number FROM patient")
all_register <- paste0(all_patients$register_number, collapse = ",")


all_person <- dbGetQuery(con, sprintf(
  "select register_number, h.name as hospital, li.code as gender, nip 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  where register_number = any(array[%s])", all_register)) %>%
  mutate(SEX = ifelse(gender == "SEXO_HOMBRE", "M", "F")) %>%
  as_tibble()



hemogram <- dbGetQuery(con, 
                       sprintf("select register_number, blasts_percent, leukocytes,
        neutrophils, hemoglobin, platelets, monos_percent  
        from patient p 
        join analytical_data ad on p.id = ad.patient_id 
        join hemogram_data hd on ad.id = hd.hemogram_data_id 
        where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
  as_tibble() %>%
  mutate(MONOCYTES = leukocytes*monos_percent/100,
         PB_BLAST = blasts_percent,
         WBC = leukocytes,
         ANC = neutrophils,
         HB = hemoglobin,
         PLT = platelets)  %>%
  select(register_number, MONOCYTES, PB_BLAST, WBC, ANC, HB, PLT)


morphological <- dbGetQuery(con, sprintf(
  "select register_number, blastos_mo
  from patient p 
  join analytical_data ad on p.id = ad.patient_id 
  join morphological_data md on ad.id = md.morphologycal_data_id 
  where register_number in (%s) and analysis_time = 1045", all_register)) %>%
  as_tibble() %>%
  mutate(BM_BLAST = blastos_mo) %>%
  select(register_number, BM_BLAST)

## Remove a copy of an individual
morphological <- morphological[!duplicated(morphological$register_number), ]

## Biochemistry
biochem_data <-   dbGetQuery(con, sprintf(
  "select register_number, erythropoietin     
  from patient p 
  join analytical_data ad on p.id = ad.patient_id 
  join biochemistry_data b on ad.id = b.bio_data_id 
  where register_number in (%s) and analysis_time = 1045"
  , all_register)) %>%
  as_tibble()


## Pronostico ####

all_pronostico <- dbGetQuery(con, sprintf(
  "select id, register_number, ipss_mol_risk_group, ipssr_risk_group, ipss, ipssr,
  ipss_mol, ipss_calculation_status 
  from patient 
  join prognostic_scoring ps on patient.id = ps.diagnosis_data 
  where register_number = any(array[%s])", all_register)) %>%
  as_tibble()

ipssm <- unique(all_pronostico$ipss_mol_risk_group  )
ipssm <- ipssm[!is.na(ipssm)]
codigos_ipssm <- dbGetQuery(con, 
                            paste("select id, code from list_item where id in (", 
                                  paste(ipssm, collapse = ","), ")")
) %>%
  mutate(IPSSM = code,
         IPSSM =  recode(IPSSM, "IPSS_MOL_GRUPO_RIESGO_MUY_BAJO" = "Very-Low",
                         "IPSS_MOL_GRUPO_RIESGO_BAJO" = "Low",
                         "IPSS_MOL_GRUPO_RIESGO_MOD_BAJO" = "Moderate-Low",
                         "IPSS_MOL_GRUPO_RIESGO_MOD_ALTO" = "Moderate-High",
                         "IPSS_MOL_GRUPO_RIESGO_ALTO" = "High",
                         "IPSS_MOL_GRUPO_RIESGO_MUY_ALTO" = "Very-High")) %>%
  select(-code) 


ipssr <- unique(all_pronostico$ipssr_risk_group  )
ipssr <- ipssr[!is.na(ipssr)]
codigos_ipssr <- dbGetQuery(con, 
                            paste("select id, code from list_item where id in (", 
                                  paste(ipssr, collapse = ","), ")")
) %>%
  mutate(IPSSR = code,
         IPSSR =  recode(IPSSR, "IPSSR_GRUPO_RIESGO_MUY_BAJO" = "Very-Low",
                         "IPSSR_GRUPO_RIESGO_BAJO" = "Low",
                         "IPSSR_GRUPO_RIESGO_INT" = "Int",
                         "IPSSR_GRUPO_RIESGO_ALTO" = "High",
                         "IPSSR_GRUPO_RIESGO_MUY_ALTO" = "Very-High")) %>%
  select(-code) 

all_pronostico2 <- left_join(all_pronostico, codigos_ipssm, 
                             by = join_by(ipss_mol_risk_group == id)) %>%
  left_join(codigos_ipssr, by = join_by(ipssr_risk_group == id))




mutations_sum <- dbGetQuery(con,
                            sprintf("select register_number, has_mutations, sample_type, ngs_done,
    ngs_board, name as ngs_name, genes as ngs_genes
    from patient p 
    join analytical_data ad on p.id = ad.patient_id 
    join molecular_data md on md.molec_data_id = ad.id
    join ngs_board nb on nb.id = md.ngs_board
    where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
  as_tibble() %>%
  mutate(has_mutations = ifelse(has_mutations == 36094, "YES", 
                                ifelse(has_mutations == 36095, "NO", has_mutations)),
         ngs_done = ifelse(ngs_done == 36094, "YES", 
                           ifelse(ngs_done == 36095, "NO", ngs_done)))

mutations <- dbGetQuery(con,
                        sprintf("select register_number, gen, type, vaf, status
    from patient p 
    join analytical_data ad on p.id = ad.patient_id 
    join mutation m on m.molecular_data_id = ad.id
    where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
  as_tibble() 

mutations_epe <- filter(mutations, register_number %in% treatment$register_number)
treatment <- dbGetQuery(con, "SELECT p.register_number, tdd.drug_name_id, tdd.treatment_data_id, tdd.*, td.*
  FROM treatment_drug_data tdd
LEFT JOIN treatment_data td ON tdd.treatment_data_id = td.id
RIGHT JOIN patient p ON td.patient_id = p.id
WHERE tdd.drug_name_id IN (36090, 36089, 36085, 36053, 36032, 4607613)")


treatment_present <- mutations_sum %>% filter(register_number %in% treatment$register_number)


epe_data2 <- left_join(all_person, hemogram, by = "register_number") %>%
  left_join(morphological, by = "register_number") %>%
  left_join(all_pronostico2, by ="register_number") %>%
  left_join(biochem_data, by = "register_number") %>%
  filter(register_number %in% treatment$register_number)

save(epe_data2, treatment, mutations_epe, file = "epe_data2_GESMD_v2.Rdata")
