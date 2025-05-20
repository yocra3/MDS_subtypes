library(tidyverse)
library(lubridate)
#library(ssh)
library(sys)
library(RPostgres)
library(ipssm)

#primero conectar tunel SSH via fichero bash
ssh -i private.pem -N -L 5556:prehm-pro-rds-i.cfeequry4htk.eu-west-3.rds.amazonaws.com:5432 ec2-user@35.181.32.109

drv <- RPostgreSQL::PostgreSQL()
con <- dbConnect(RPostgres::Postgres(), user="cruizarenas_resmdro",
                 password='s"4/Z7bG=m%9j3W',
                 dbname="masihdas_resmd",
                 port=5556, # el puerto es el mismo al que rediriges el tunel ssh
                 host="localhost")
#prueba para ver que bases de datos hay y que funciona el tunel
tabs <- dbListTables(con)
fields <- lapply(tabs, dbListFields, conn = con)
fields_df <- Reduce(rbind, lapply(fields, function(x) tibble(field = x))) %>%
  mutate(table = rep(tabs, lengths(fields)))



## Get patients ####
all_patients <- dbGetQuery(con, "SELECT id, register_number FROM patient")


all_register <- paste0(all_patients$register_number, collapse = ",")

all_person <- dbGetQuery(con, sprintf(
  "select register_number, h.name as hospital, li.code as gender, birthdate 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  where register_number = any(array[%s])", all_register)) %>%
  mutate(SEX = ifelse(gender == "SEXO_HOMBRE", "M", "F")) %>%
  as_tibble()

all_diagnosis <- dbGetQuery(con, sprintf(
  "select p.id, register_number, birthdate, diagnosis_date, li1.code as smd_secundario,
  li2.code as who2017, li3.code as who2008
  from patient p 
  join diagnosis_data dd on p.id = dd.id
  left join list_item li1 on dd.smd_id = li1.id
  left join list_item li2 on who2017_id = li2.id
  left join list_item li3 on who2008_id = li3.id
  where register_number = any(array[%s])", all_register)) %>%
  as_tibble() %>% 
  mutate(AGE = time_length(difftime(diagnosis_date, birthdate), "years"))

all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_SINDROME_5Q"] <- "WHO2017_SMD_DEL_5Q"
all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_LMMC"] <- "WHO2017_LMMCX"
all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_CRDM"] <- "WHO2017_SMD_DM"
all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_AREB_1"] <- "WHO2017_SMD_EB_1"
all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_AREB_2"] <- "WHO2017_SMD_EB_2"
all_diagnosis$who2017[is.na(all_diagnosis$who2017) & all_diagnosis$who2008=="WHO2008_ARS"] <- "WHO2017_SMD_SA_DU"


## Estado ####
estado_actual <- dbGetQuery(con, 
  sprintf("select register_number, assessment_date, death_reason,
  status, disease_state_death 
  from patient pr 
  join current_status cs on pr.id = cs.id 
  where register_number = any(array[%s])", all_register))  %>%
  as_tibble() 

status <- unique(estado_actual$status)
status <- status[!is.na(status)]
codigos_status <- dbGetQuery(con, 
  paste("select id, code from list_item where id in (", 
    paste(status, collapse = ","), ")")
) %>%
  mutate(estado = code) %>%
  select(-code)

muerte <- unique(estado_actual$death_reason  )
muerte <- muerte[!is.na(muerte)]
codigos_muerte <- dbGetQuery(con, 
                             paste("select id, code from list_item where id in (", 
                                   paste(muerte, collapse = ","), ")")
) %>%
  mutate(causa_muerte = code) %>%
  select(-code)


estado_actual2 <- left_join(estado_actual, codigos_status, 
                            by = join_by(status == id)) %>%
  left_join(codigos_muerte, by = join_by(death_reason == id))

### Hay pacientes con causa de muerte sin estado muerto ¿?
patient_basic <- left_join(all_person, all_diagnosis, by = "register_number" ) %>%
  left_join(estado_actual2, by = "register_number" ) %>%
  mutate(OS_YEARS = time_length(difftime(assessment_date , diagnosis_date), "years"),
         OS_STATUS = ifelse(estado == "ESTADO_PACIENTE_MUERTO", 1, 
                            ifelse(is.na(estado), NA, 0)))


## Pronostico ####

all_pronostico <- dbGetQuery(con, sprintf(
  "select id, register_number, ipss_mol_risk_group, ipssr_risk_group, ipss, ipssr, ipss_calculation_status 
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
  mutate(IPSSM = code) %>%
  select(-code)

ipssr <- unique(all_pronostico$ipssr_risk_group  )
ipssr <- ipssr[!is.na(ipssr)]
codigos_ipssr <- dbGetQuery(con, 
                            paste("select id, code from list_item where id in (", 
                                  paste(ipssr, collapse = ","), ")")
) %>%
  mutate(IPSSR = code) %>%
  select(-code)

all_pronostico2 <- left_join(all_pronostico, codigos_ipssm, 
                            by = join_by(ipss_mol_risk_group == id)) %>%
  left_join(codigos_ipssr, by = join_by(ipssr_risk_group == id))

## Muy pocas muestras con IPSS-mol

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


### Karyotype ####
cariotype_raw <- dbGetQuery(con,
  sprintf("select register_number, cariotype_description, ipssr_cytogenetic 
    from patient p 
    join analytical_data ad on p.id = ad.patient_id 
    join citogenetical_data cd on cd.cito_data_id = ad.id 
    where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
    as_tibble()


ipssr_cyto <- unique(cariotype_raw$ipssr_cytogenetic  )
ipssr_cyto <- ipssr_cyto[!is.na(ipssr_cyto)]
codigos_ipssr_cyto <- dbGetQuery(con, 
                            paste("select id, code from list_item where id in (", 
                                  paste(ipssr_cyto, collapse = ","), ")")
) %>%
  mutate(CYTO_IPSSR = code,
         CYTO_IPSSR =  recode(CYTO_IPSSR, "IPSSR_CITOGENETICO_MUY_POBRE" = "Very Poor",
                              "IPSSR_CITOGENETICO_POBRE" = "Poor",
                              "IPSSR_CITOGENETICO_INTERMEDIO" = "Intermediate",
                              "IPSSR_CITOGENETICO_BUENO" = "Good",
                              "IPSSR_CITOGENETICO_MUY_BUENO" = "Very Good")) %>%
  select(-code) 
 
  
cariotype_raw2 <- left_join(cariotype_raw, codigos_ipssr_cyto, 
                             by = join_by(ipssr_cytogenetic == id)) 

cariotype <- dbGetQuery(con,
  sprintf("select register_number, code, has_anomaly, cito_data_id 
    from patient p 
    join analytical_data ad on p.id = ad.patient_id 
    join citogenetical_data_anomalies_normalized cd on cd.cito_data_id = ad.id 
    where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
  as_tibble()

top_car_df <- cariotype %>%
  group_by(code) %>%
  summarize(N = sum(has_anomaly)) %>%
  arrange(desc(N))

top_car <- top_car_df$code[1:10]
names(top_car) = top_car

lapply(top_car, function(car){
  filter(cariotype_raw, 
         register_number %in% filter(cariotype, code == car & has_anomaly)$register_number)
  
})

mapping <- tibble(code = paste0("ANOMALIAS_CITOGENETICAS_", c("XXXIV", "IV", "II", "XXX", "X", "XXIX", "XXXV", "VI")),
                  event = c("del5q", "plus8", "delY", "del7", "del20q", "complex", "del7q", "del11q"))


cariotype_tab <- cariotype %>%
  right_join(mapping, by = "code") %>%
  mutate(value = ifelse(has_anomaly, 1, 0)) %>%
  select(register_number, event, value) %>%
  pivot_wider(names_from = event, values_from = value) %>%
  mutate(del7_7q = ifelse(del7 == 1 | del7q == 1, 1, 0),
         del17_17p = NA) %>% 
  left_join(select(cariotype_raw2, register_number, CYTO_IPSSR), by = "register_number" )

# Mutations ####

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

## Descartamos las que no tienen has_mutations
mutations <- dbGetQuery(con,
                        sprintf("select register_number, gen, type, vaf, status
    from patient p 
    join analytical_data ad on p.id = ad.patient_id 
    join mutation m on m.molecular_data_id = ad.id
    where analysis_time = 1045 and register_number in (%s)", all_register)) %>%
  as_tibble() 

muts <- unique(mutations$status  )
muts <- muts[!is.na(muts)]
codigos_muts <- dbGetQuery(con, 
                           paste("select id, code from list_item where id in (", 
                                 paste(muts, collapse = ","), ")")
) %>%
  mutate(mutation = code) %>%
  select(-code)


mutations2 <- left_join(mutations, codigos_muts, 
                        by = join_by(status == id)) %>%
  mutate(value = ifelse(mutation == "MUTATION_STATUS_NO_MUTADO", 0, 1))

mutations_tab <- select(mutations2, register_number, gen, value) %>%
  group_by(register_number, gen) %>%
  summarize(mut = sum(value)) %>%
  mutate(mut = ifelse(!gen %in% c("TP53", "TET2") & mut > 1, 1, mut)) %>%
  pivot_wider(names_from = gen, values_from = mut) %>%
  mutate(TET2 = ifelse(TET2 == 0, 0, 1), 
         TP53 = ifelse(TP53 == 0, 0, 1),
         TP53mut = TP53, 
         TP53loh = NA,
         MLL_PTD = NA) %>%
  left_join(mutations2 %>%
              filter(gen %in% c("TP53", "TET2")) %>%
              group_by(register_number, gen) %>%
              mutate(vaf = ifelse(is.na(vaf), 0, vaf)) %>%
              summarize(vaf = sum(vaf, na.rm = TRUE)/100) %>%
              pivot_wider(names_from = gen, values_from = vaf) %>%
              mutate(TET2bi = ifelse(TET2 >= 0.5, 1, 0), 
                     TET2other = ifelse(TET2 > 0 & TET2 < 0.5, 1, 0),
                     TP53multi = ifelse(TP53 >= 0.5, 1, 0)) %>%
              select(-TET2, -TP53),
            by = "register_number")



## Select mutations for this study
sel_muts <- c("TET2", "ASXL1", "SRSF2", "DNMT3A", "RUNX1", 
              "STAG2", "U2AF1", "EZH2", "ZRSR2", "TET2bi", "TET2other",
              "SF3B1", "TP53multi", "FLT3", "BCOR", "BCORL1", "CBL",
              "CEBPA",	"ETV6",	"IDH1",	"IDH2", "KRAS", "NF1", "NPM1",
              "NRAS",	"SETBP1",	"ETNK1", "GATA2",	"GNB1",	"PHF6",	"PPM1D", "PRPF8",
              "PTPN11",	"WT1", "MLL_PTD")
mutations_sel <- select(mutations_tab, register_number, all_of(sel_muts), starts_with("TP53"))

  
mutations2 %>% 
  group_by(register_number) %>% 
  summarize(N_mut = sum(mutation == "MUTATION_STATUS_MUTADO")) %>%
  arrange(N_mut)

## Final merge
gesmd_data <- left_join(patient_basic, hemogram, by = "register_number") %>%
  left_join(morphological, by = "register_number") %>%
  left_join(all_pronostico2, by = "register_number") %>%
  left_join(cariotype_tab, by = "register_number") %>%
  left_join(mutations_sel, by = "register_number") %>%
  mutate(consensus = ifelse(TP53multi == 1 & BM_BLAST <= 20, "Mutated TP53",
                   ifelse(del5q == 1 & del7q == 0 & BM_BLAST <= 5, "del5q",
                          ifelse(SF3B1 > 0 & del7q == 0 & complex == 0 & BM_BLAST <= 5, "mutated SF3B1",
                                 ifelse(BM_BLAST <= 5, "Low blasts",
                                        ifelse(BM_BLAST > 10, "MDS-IB2",
                                               ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other")))))))
save(gesmd_data, file = "results/gesmd_data_all.Rdata")

write.csv(gesmd_data, file = "results/gesmd_data.csv")

## Compute IPSSM ####
ipssm_raw <- IPSSMread("results/gesmd_data.csv")
ipssm_process <- IPSSMprocess(ipssm_raw)
ipssm_res <- IPSSMmain(ipssm_process)
ipssm_annot <- IPSSMannotate(ipssm_res)
## No funciona. Falta mucha información

gesmd_low <- filter(gesmd_data, consensus == "Low blasts")
save(gesmd_low, file = "results/gesmd_data_low.Rdata")

## Match with Raul
raul_data <- read_csv("pacientesngsvcf.csv") %>%
  mutate(selected = ifelse(patient_id %in% gesmd_data$register_number, "Included", "Missing")) %>%
  left_join(select(gesmd_data, register_number, consensus) %>%
              mutate(register_number = as.double(register_number)),
            by = join_by(patient_id == register_number))

table(raul_data$ngs_done, raul_data$selected)
table(raul_data$ngs_done, !is.na(raul_data$consensus))
