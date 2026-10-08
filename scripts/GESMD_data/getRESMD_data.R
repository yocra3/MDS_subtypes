library(tidyverse)
library(lubridate)
library(readxl)
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

## Evolution ####
evolution <- dbGetQuery(con, 
                            sprintf("select register_number, who2017_type_id, evol_date
  from patient pr 
  join evolution_data ed on pr.id = ed.id 
  where register_number = any(array[%s])", all_register))  %>%
  as_tibble() 

evolution_type <- unique(evolution$who2017_type_id)
evolution_type <- evolution_type[!is.na(evolution_type)]
codigos_evo <- dbGetQuery(con, 
                             paste("select id, code from list_item where id in (", 
                                   paste(evolution_type, collapse = ","), ")")
) %>%
  mutate(who2017_evo = code) %>%
  select(-code)

evolution2 <- left_join(evolution, codigos_evo, 
                            by = join_by(who2017_type_id == id))

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
  left_join(evolution2, by = "register_number") %>%
  mutate(OS_YEARS = time_length(difftime(assessment_date , diagnosis_date), "years"),
         OS_STATUS = ifelse(estado == "ESTADO_PACIENTE_MUERTO", 1, 
                            ifelse(is.na(estado), NA, 0)),
         AMLt_STATUS = ifelse(is.na(who2017_evo) | who2017_evo != "WHO2017_LMA", 0, 1),
         evol_date2 = ifelse(is.na(evol_date), assessment_date, evol_date),
         evol_date2 = as.Date(evol_date2),
         AMLt_YEARS = time_length(difftime(evol_date2 , diagnosis_date), "years"))


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

## Muy pocas muestras con IPSS-mol
## Hemogram ####
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


## Karyotype ####

## Read mapping
kar_map <- read_xlsx("data/GESMD_genetics/anomalias citogenéticas Carlos.xlsx")

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

# top_car_df <- cariotype %>%
#   group_by(code) %>%
#   summarize(N = sum(has_anomaly)) %>%
#   arrange(desc(N))
# 
# top_car <- top_car_df$code[1:10]
# names(top_car) = top_car
# 
# lapply(top_car, function(car){
#   filter(cariotype_raw, 
#          register_number %in% filter(cariotype, code == car & has_anomaly)$register_number)
#   
# })
# 
# mapping <- tibble(code = paste0("ANOMALIAS_CITOGENETICAS_", c("XXXIV", "IV", "II", "XXX", "X", "XXIX", "XXXV", "VI")),
#                   event = c("del5q", "plus8", "delY", "del7", "del20q", "complex", "del7q", "del11q"))


cariotype_tab <- cariotype %>%
  right_join(kar_map, by = join_by(code == variables)) %>%
  mutate(value = ifelse(has_anomaly, 1, 0),
         event = case_when(
           anomalias == "+8" ~ "plus8",
           anomalias == "-7" ~ "del7",
           anomalias == "-Y" ~ "delY",
           anomalias == "-5" ~ "del5",
           anomalias == "-13" ~ "del13",
           anomalias == "+19" ~ "plus19",
           TRUE ~ str_replace_all(anomalias, "[()]", "")
         )) %>%
  select(register_number, event, value) %>%
  pivot_wider(names_from = event, values_from = value) %>%
  mutate(N_aberrations = rowSums(pick(-1)),
    complex = ifelse(N_aberrations >= 3, 1, 0),
        del7_7q = ifelse(del7 == 1 | del7q == 1, 1, 0),
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

tp53_tet2_tab <- mutations2 %>%
  filter(gen %in% c("TP53", "TET2")) %>%
  group_by(register_number, gen) %>%
  mutate(vaf = ifelse(is.na(vaf), 0, vaf)) %>%
  summarize(n_rel = sum(vaf >= 10), 
            vaf = sum(vaf, na.rm = TRUE)/100) %>%
  group_by(register_number) 

tet2_tab <- subset(tp53_tet2_tab, gen == "TET2") %>%
  mutate(TET2bi = ifelse(vaf >= 0.5, 1, 0), 
         TET2other = ifelse(vaf > 0 & vaf < 0.5, 1, 0))
tp53_tab <- subset(tp53_tet2_tab, gen == "TP53") %>%
         mutate(TP53multi = ifelse(vaf >= 0.5 | n_rel >= 2 , 1, 0),
                TP53mono = ifelse(vaf > 0 & TP53multi == 0, 1, 0))
  

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
  left_join(select(tet2_tab, register_number, TET2bi, TET2other), by = "register_number") %>%
  left_join(select(tp53_tab, register_number, TP53multi, TP53mono), by = "register_number")


## Check TET2
tet2 <- mutations2 %>%
  filter(gen == "TET2") %>%
  group_by(register_number, gen) %>%
  mutate(vaf = ifelse(is.na(vaf), 0, vaf)) %>%
  summarize(comp_vaf = sum(vaf, na.rm = TRUE)/100,
            n_mut = sum(vaf > 0),
            max_vaf = max(vaf/100, na.rm = TRUE)) %>%
  filter(n_mut > 0) %>%
  ungroup() %>%
  mutate(class = case_when(
    n_mut == 1 & max_vaf < 0.5 ~ "TET2mono",
    n_mut == 1 & max_vaf > 0.5 ~ "TET2bi",
    n_mut > 1 & comp_vaf < 0.5 ~ "TET2mono",
    n_mut > 1 & max_vaf > 0.5 ~ "TET2bi",
    n_mut > 1 & comp_vaf > 0.5 & max_vaf < 0.5 ~ "TET2bi - doubts"
  ))

tet2 %>%
  filter(n_mut > 1) %>%
  ggplot(aes(x = comp_vaf, y = max_vaf, col = factor(n_mut))) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5)

summarize(tet2, 
          all = n(),
          clear = sum(n_mut == 1),
          tet2bi = sum(n_mut > 1 & max_vaf > 0.5),
          tet2m = sum(n_mut > 1 & comp_vaf < 0.5),
          unknown = sum(n_mut > 1 & comp_vaf > 0.5 & max_vaf < 0.5),
          p_unknown = mean(n_mut > 1 & comp_vaf > 0.5 & max_vaf < 0.5))


filter(tet2, n_mut > 1 & comp_vaf > 0.5 & max_vaf < 0.5)



## Select mutations for this study
sel_muts <- c("TET2", "ASXL1", "SRSF2", "DNMT3A", "RUNX1", 
              "STAG2", "U2AF1", "EZH2", "ZRSR2", "TET2bi", "TET2other",
              "SF3B1","TP53mono","TP53multi", "FLT3", "BCOR", "BCORL1", "CBL",
              "CEBPA",	"ETV6",	"IDH1",	"IDH2", "KRAS", "NF1", "NPM1",
              "NRAS",	"SETBP1",	"ETNK1", "GATA2",	"GNB1",	"PHF6",	"PPM1D", "PRPF8",
              "PTPN11",	"WT1", "MLL_PTD")
mutations_sel <- select(mutations_tab, register_number, all_of(sel_muts), starts_with("TP53"))

  
mutations2 %>% 
  group_by(register_number) %>% 
  summarize(N_mut = sum(mutation == "MUTATION_STATUS_MUTADO")) %>%
  arrange(N_mut)


# Treatment ####
trasplante <- dbGetQuery(con, sprintf(
  "select register_number, transplanted from patient 
  join tph_data on patient.id = tph_data.id where register_number
  in (%s) order by register_number", all_register)) 

transfusion <- dbGetQuery(con, sprintf(
  "select register_number from patient p join 
  transfusion t on p.id = t.patient_id where register_number in (%s)
  order by register_number", all_register))  %>%
  distinct() %>%
  mutate(transfusion = TRUE)

tratamiento<-dbGetQuery(con, sprintf(
  "select register_number, drug_name_id, support_drug_name_id,
   start_date, end_date, treatment_line, description
   from patient p join treatment_data td on p.id = td.patient_id
   left join treatment_drug_data tdd on td.id = tdd.treatment_data_id where
   register_number in (%s) order by register_number", 
  all_register)) 


tratamiento$register_number<-as.numeric(tratamiento$register_number)

tratamiento2 <- dbGetQuery(con, 
                         sprintf("select register_number, drug_name_id, support_drug_name_id,
          td.start_date, end_date, treatment_line, description
          from patient p join treatment_data td on p.id = td.patient_id
          join treatment_cycle_data tcd on td.id = tcd.treatment_data_id left join
          treatment_drug_data tdd2 on tcd.id = tdd2.cycle_data_id where
          register_number in (%s) order by register_number", 
                                 all_register)) 

tratamiento2$register_number<-as.numeric(tratamiento2$register_number)

t3 <- rbind(tratamiento, tratamiento2) %>% group_by(register_number) %>% 
  arrange(start_date, .by_group = T)

t4 <- unite(t3, drug_id, c(drug_name_id, support_drug_name_id))
t4$drug_id <- gsub("NA_","", t4$drug_id)
t4$drug_id <- gsub("_NA", "", t4$drug_id)
t4$drug_id[t4$drug_id == "NA"] <- NA
t4$drug_id <- parse_integer(t4$drug_id)

lista_tratamientos <- dbGetQuery(con, sprintf(
  "select id, code as tratamiento from list_item where id in (%s) order by id",
  paste0(unique(t4$drug_id[!is.na(t4$drug_id)]), collapse = ",")))

lista_tratamientos$id <- as.numeric(lista_tratamientos$id)

treatment_merge <- left_join(t4, lista_tratamientos, by = join_by(drug_id == id)) %>%
  mutate(hma = ifelse(tratamiento %in% c("NOM_FARMACO_AZACITIDINA", "NOMBRE_FARMACO_DECITABINA", "NOM_FARMACO_AGENTES_HIPOMETILANTES"), 1, 0),
         azacitidina = ifelse(tratamiento == "NOM_FARMACO_AZACITIDINA", 1, 0),
         lenalidomida = ifelse(tratamiento == "NOM_FARMACO_LENALIDOMIDA", 1, 0)) %>%
  select(register_number, hma, azacitidina, lenalidomida) %>%
  pivot_longer(cols = -register_number) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  group_by(register_number, name) %>%
  summarize(value = max(value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(register_number = as.integer(register_number))
           

## Final merge ####
gesmd_data <- left_join(patient_basic, hemogram, by = "register_number") %>%
  left_join(morphological, by = "register_number") %>%
  left_join(all_pronostico2, by = "register_number") %>%
  left_join(cariotype_tab, by = "register_number") %>%
  left_join(mutations_tab, by = "register_number") %>%
  left_join(treatment_merge, by = "register_number") %>%
  left_join(trasplante, by = "register_number") %>%
  left_join(transfusion, by = "register_number") %>%
  mutate(consensus = ifelse(if_any(c(TP53, del5q, SF3B1, BM_BLAST), is.na), "Undetermined",
    ifelse((TP53multi == 1 | (TP53mono == 1 & del17p13 == 1)) & BM_BLAST <= 20, "MDS-TP53",
           ifelse(del5q == 1 & del7q == 0 & del7 == 0 & BM_BLAST <= 5 & N_aberrations < 3 & complex == 0, "MDS-del5q",
                  ifelse(SF3B1 == 1 & del7q == 0 & del7 == 0 & RUNX1 == 0 & complex == 0 & BM_BLAST <= 5, "MDS-SF3B1",
                         ifelse(BM_BLAST <= 5, "MDS-LB",
                                ifelse(BM_BLAST > 10, "MDS-IB2",
                                       ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other"))))))
  ))

gesmd_data$register_number <- as.character(gesmd_data$register_number)
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

## Código Irene ####
all_person_irene <- dbGetQuery(con, sprintf(
  "select register_number, nip, h.name as hospital, li.code as gender, birthdate 
  from patient p 
  join hospital h on p.hospital_id = h.id 
  join list_item li on p.gender_id = li.id 
  where register_number = any(array[%s])", all_register)) %>%
  mutate(SEX = ifelse(gender == "SEXO_HOMBRE", "M", "F")) %>%
  as_tibble()



cun <- subset(all_person_irene, hospital == "Clinica Universidad de Navarra") %>%
  left_join(gesmd_data, by = "register_number")

library(readxl)
library(writexl)

data1 <- read_xlsx("./data/Pacientes_Irene.xlsx", sheet = "Hoja1") %>%
  mutate(nip = as.character(`Nº de historia`))

data1_merge <- left_join(data1, cun, by = "nip")

data2 <- read_xlsx("./data/Pacientes_Irene.xlsx", sheet = "Solo tto CUN o no tto") %>%
  mutate(nip = as.character(`Nº de historia`))

data2_merge <- left_join(data2, cun, by = "nip")

write_xlsx(list(Hoja1 = data1_merge, "Solo tto CUN o no tto" = data2_merge), 
           "./data/Pacientes_Irene_merge.xlsx")
