library(tidyverse)
library(lubridate)
#library(ssh)
library(sys)
library(RPostgres)
library(writexl)
library(dplyr)
library(RCytoGPS)

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



## Get patients
all_patients <- dbGetQuery(con, "SELECT register_number FROM patient")


all_register <- all_patients$register_number

all_person <- dbGetQuery(con, sprintf(
  "select register_number, h.name as hospital, 
  li.code as gender, birthdate from patient p join hospital h on
  p.hospital_id = h.id join list_item li on p.gender_id = li.id where 
  register_number = any(array[%s]) order by register_number",
  paste0(all_register, collapse = ","))) %>%
  mutate(sex = ifelse(gender == "SEXO_HOMBRE", "M", "F"))

all_diagnosis <- dbGetQuery(con, sprintf(
  "select register_number, birthdate, diagnosis_date, li1.code as smd_secundario
  from patient p 
  join diagnosis_data dd on p.id = dd.id
  left join list_item li1 on dd.smd_id = li1.id
  where register_number in (%s)",
  paste0(all_num_web, collapse = ","))) %>%
  as_tibble()



#### Aqui empezamos a seleccionar pacientes, editar segun necesidades####

patients <- dbGetQuery(con, 
                       "select register_number from patient where 
   register_number = any(array[8631,11777,11863,12032,13123,13321,13957,14945,14956,16069,16337,16353,16356,16369,16372,16380,16387,16413,16414,16903,16905,17392,17418,17422,17428,17433,17444,17584,17627,17846,17868,17872,17873,17874,18076,18088,18216,18262,18266,18278,18331,18378,18379,18385,18427,18448,18449,18557,18580,18652,18654,18655,18656,18659,18661,18662,18663,18664,18667,18668,18669,18670,18671,18674,18675,18676,18677,18679,18680,18682,18683,18685,18689,18691,18692,18693,18714,18748,18755,18841,18874,18878,18924,18936,18964,18970,18990,19016,19037,19040,19041,19063,19071,19132,19138,19227,15084,15444,16410,17266,17414,17416,17420,17442,17492,17493,17900,18037,18042,18139,18578,18605,18611,18612,18613,18615,18660,18673,18681,18686,18688,18764,18926,19039,19048,19066,19081,19093,19128,19329,19444,19466])"
)


all_num_web <- unique(as.numeric(patients$register_number))
df_num_web <- data.frame(num_web = all_num_web)


register_status <- dbGetQuery(con, sprintf(
  "select register_number, useful_id from patient p join
  register_status cs on p.id = cs.id where register_number in (%s)",
  paste0(all_num_web, collapse = ",")))

table(register_status$useful_id, useNA = "always")

#### patient: id, centro, country, gender, birth date ####

person <- dbGetQuery(con, sprintf(
  "select register_number, h.name as hospital, 'Spain' as country,
  li.code as gender, birthdate from patient p join hospital h on
  p.hospital_id = h.id join list_item li on p.gender_id = li.id where 
  register_number = any(array[%s]) order by register_number",
  paste0(all_num_web, collapse = ",")))


nhc <- dbGetQuery(con, sprintf(
  "select register_number, nip, identification from patient p where 
  register_number in (%s) order by register_number",
  paste0(all_num_web, collapse = ",")))

write.table(nhc, file = "num_web_y_nhc.csv", sep = ";",
            row.names = F)

person$gender[person$gender == "SEXO_HOMBRE"] <- "M"
person$gender[person$gender == "SEXO_MUJER"] <- "F"

person$register_number <- as.numeric(person$register_number)

#### diagnosis_data: date, who2017, age at diagnosis ####

diagnosis <- dbGetQuery(con, sprintf(
  "select register_number, birthdate, diagnosis_date, li3.code as who2008, li1.code as who2017,
  li2.code as smd_secundario
  from patient p join diagnosis_data dd on p.id = dd.id left join list_item li1 on
  who2017_id = li1.id left join list_item li2 on dd.smd_id = li2.id
  left join list_item li3 on who2008_id =li3.id where 
  register_number in (%s) order by register_number",
  paste0(all_num_web, collapse = ",")))


#solo ponemos el diagnostico 2017 usando el 2008 si esta vacio
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_SINDROME_5Q"] <- "WHO2017_SMD_DEL_5Q"
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_LMMC"] <- "WHO2017_LMMCX"
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_CRDM"] <- "WHO2017_SMD_DM"
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_AREB_1"] <- "WHO2017_SMD_EB_1"
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_AREB_2"] <- "WHO2017_SMD_EB_2"
diagnosis$who2017[is.na(diagnosis$who2017) & diagnosis$who2008=="WHO2008_ARS"] <- "WHO2017_SMD_SA_DU"

diagnosis$who2022[diagnosis$who2017=="WHO2017_LMMC"] <- "WHO2022_LMMCX"
diagnosis$who2017[diagnosis$who2017=="WHO2017_LMMC"] <- "WHO2017_LMMCX"


diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_DEL_5Q"] <- "WHO2022_SMD_DEL_5Q"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_DM"] <- "WHO2022_SMD_LB"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_DU"] <- "WHO2022_SMD_LB"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_SA_DU"] <- "WHO2022_SMD_SF3B1" 
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_SA_DM"] <- "WHO2022_SMD_SF3B1"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_EB1"] <- "WHO2022_SMD_IB_1"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_EB_2"] <- "WHO2022_SMD_IB_2" 
diagnosis$who2022[diagnosis$who2017=="WHO2017_LMMCX"] <- "WHO2022_LMMCX" 
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_SMP_SA_T"] <- "WHO2022_SMD_SMP_SA_SF3B1"
diagnosis$who2022[diagnosis$who2017=="WHO2017_LMMC_0"] <- "WHO2022_LMMC_1"
diagnosis$who2022[diagnosis$who2017=="WHO2017_LMMC_1"] <- "WHO2022_LMMC_1"
diagnosis$who2022[diagnosis$who2017=="WHO2017_LMMC_2"] <- "WHO2022_LMMC_2"
diagnosis$who2022[diagnosis$who2017=="WHO2017_LMA"] <- "WHO2022_LMA"
diagnosis$who2022[diagnosis$who2017=="WHO2017_SMD_SMP_INCLASIFICABLE"] <- "WHO2022_SMD_SMP_INCLASIFICABLE"

diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_DEL_5Q"] <- "ICC2022_SMD_DEL_5Q"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_DM"] <- "ICC2022_SMD_NOS_DM"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_DU"] <- "ICC2022_SMD_NOS_DU"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_SA_DU"] <- "ICC2022_SMD_SF3B1" 
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_SA_DM"] <- "ICC2022_SMD_SF3B1"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_EB1"] <- "ICC2022_SMD_EB"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_EB_2"] <- "ICC2022_SMD_LMA" 
diagnosis$icc2022[diagnosis$who2017=="WHO2017_LMMCX"] <- "ICC2022_LMMCX" 
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_SMP_SA_T"] <- "ICC2022_SMD_SMP_SA_SF3B1"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_LMMC_0"] <- "ICC2022_LMMC_1"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_LMMC_1"] <- "ICC2022_LMMC_1"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_LMMC_2"] <- "ICC2022_LMMC_2"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_LMA"] <- "ICC2022_LMA"
diagnosis$icc2022[diagnosis$who2017=="WHO2017_SMD_SMP_INCLASIFICABLE"] <- "ICC2022_SMD_SMP_INCLASIFICABLE"


diagnosis <- diagnosis %>% mutate(edad_diagnostico = 
                                    floor(abs(time_length(difftime(birthdate, 
                                                                   diagnosis_date), "years"))))


#diagnosis$who2017 <- gsub("WHO2017_", "", diagnosis$who2017)
#diagnosis$who2022 <- gsub("WHO2022_", "", diagnosis$who2022)


diagnosis$register_number <- as.numeric(diagnosis$register_number)

#### pronostico ####

pronostico <- dbGetQuery(con, sprintf(
  "select register_number, ipss, ipssr, ipss_calculation_status from patient join prognostic_scoring ps on patient.id = ps.diagnosis_data where 
  register_number in (%s) order by register_number",
  paste0(all_num_web, collapse = ",")))
pronostico$register_number <- as.numeric(pronostico$register_number)
final_pronostico <- merge(df_num_web,pronostico, by.x="num_web", by.y="register_number", all.x=T,all.y=F)

#### hemogram_data: analysis_date, leuko, neutro, mono, hg, platelets ####

hemogram <- dbGetQuery(con, 
                       sprintf("select register_number, analysis_date, hemoglobin, leukocytes,
        platelets, neutrophils, monos_percent, blasts_percent
        from patient p join analytical_data ad on p.id = ad.patient_id 
        join hemogram_data hd on ad.id = hd.hemogram_data_id 
        where analysis_time = 1045 
        and register_number in (%s) ORDER BY register_number", 
                               paste0(all_num_web, collapse = ",")))

hemogram$register_number <- as.numeric(hemogram$register_number)
hemogram$monocytes <- hemogram$leukocytes*hemogram$monos_percent/100

hemogram2 <- dbGetQuery(con, 
                        sprintf("select register_number, analysis_date, hemoglobin, leukocytes,
        platelets, neutrophils, monos_percent, blasts_percent
        from patient p join analytical_data ad on p.id = ad.patient_id 
        join hemogram_data hd on ad.id = hd.hemogram_data_id 
        where register_number in (%s) ORDER BY register_number", 
                                paste0(setdiff(all_num_web,hemogram$register_number), collapse = ",")))

hemogram2$register_number <- as.numeric(hemogram2$register_number)
hemogram2$monocytes <- hemogram2$leukocytes*hemogram2$monos_percent/100

hemogram <- rbind(hemogram,hemogram2)

hemogram$HB_RS <- as.numeric(as.character(cut(hemogram$hemoglobin,breaks=c(-Inf,8,10,Inf),labels=c('1.5','1','0'),right = FALSE)))
hemogram$ANC_RS  <- as.numeric(as.character(cut(hemogram$neutrophils,breaks=c(-Inf,0.8,Inf),labels=c('0.5','0'),right = FALSE)))
hemogram$PLT_RS <- as.numeric(as.character(cut(hemogram$platelets,breaks=c(-Inf,50,100,Inf),labels=c(1,0.5,0),right = FALSE)))

hemogram$R[as.numeric(hemogram$hemoglobin) < 10 & as.numeric(hemogram$neutrophils) < 1.8] <- 0.5
hemogram$R[as.numeric(hemogram$hemoglobin) < 10 & as.numeric(hemogram$platelets) < 100] <- 0.5
hemogram$R[as.numeric(hemogram$neutrophils) < 1.8 & as.numeric(hemogram$platelets) < 100] <- 0.5
hemogram$R[as.numeric(hemogram$hemoglobin) >= 10 & as.numeric(hemogram$neutrophils) >= 1.8 & as.numeric(hemogram$platelets) < 100] <- 0
hemogram$R[as.numeric(hemogram$hemoglobin) < 10 & as.numeric(hemogram$neutrophils) >= 1.8 & as.numeric(hemogram$platelets) >= 100] <- 0
hemogram$R[as.numeric(hemogram$hemoglobin) >= 10 & as.numeric(hemogram$neutrophils) < 1.8 & as.numeric(hemogram$platelets) >= 100] <- 0
hemogram$R[as.numeric(hemogram$hemoglobin) >= 10 & as.numeric(hemogram$neutrophils) >= 1.8 & as.numeric(hemogram$platelets) >= 100] <- 0



final_hemogram <- merge(df_num_web, hemogram, by.x = "num_web",
                        by.y = "register_number", all.x = T, all.Y = F)

#### biochemistry_data: LDH, ferritin ####

biochemistry<-dbGetQuery(con, 
                         sprintf("select register_number, analysis_date, ldh, ferritin 
        from patient p join analytical_data ad on p.id = ad.patient_id 
        join biochemistry_data hd on ad.id = hd.bio_data_id 
        where analysis_time = 1045 
        and register_number in (%s) ORDER BY register_number", 
                                 paste0(all_num_web, collapse = ",")))

biochemistry$register_number<-as.numeric(biochemistry$register_number)

final_biochemistry<-merge(df_num_web, biochemistry, by.x = "num_web",
                          by.y = "register_number", all.x = T, all.Y = F)

#### morphological_data: sideroblastos_anillados, bm blasts ####

morphological<-dbGetQuery(con, sprintf(
  "select register_number, analysis_date, blastos_mo, 
  sideroblastos_anillados, mielofibrosis
  from patient p join analytical_data ad on p.id = ad.patient_id join 
  morphological_data md on ad.id = md.morphologycal_data_id where register_number
  in (%s) and analysis_time = 1045 order by register_number",
  paste0(all_num_web, collapse = ","))) %>% 
  mutate(fibrosis = as.numeric(mielofibrosis) %in% c(1110, 1115, 1120))

morphological$register_number<-as.numeric(morphological$register_number)


morphological2<-dbGetQuery(con, sprintf(
  "select register_number, analysis_date, blastos_mo, 
  sideroblastos_anillados, mielofibrosis
  from patient p join analytical_data ad on p.id = ad.patient_id join 
  morphological_data md on ad.id = md.morphologycal_data_id where register_number
  in (%s) order by register_number",
  paste0(setdiff(all_num_web,morphological$register_number), collapse = ","))) %>% 
  mutate(fibrosis = as.numeric(mielofibrosis) %in% c(1110, 1115, 1120))

morphological2$register_number<-as.numeric(morphological2$register_number)

morphological <- rbind(morphological,morphological2)

morphological$BM_BLAST_RS  <- as.numeric(as.character(cut(morphological$blastos_mo,breaks=c(-Inf,2,4.99999,10,Inf),labels=c(0,1,2,3),right = TRUE)))
morphological$BM_BLAST_S  <- as.numeric(as.character(cut(morphological$blastos_mo,breaks=c(-Inf,4,10,20,Inf),labels=c(0,0.5,1.5,2),right = TRUE)))



final_morphological<-merge(df_num_web, morphological, by.x = "num_web",
                           by.y = "register_number", all.x = T, all.Y = F)

displasias<-dbGetQuery(con, sprintf(
  "select register_number, diseritropoyesis, has_diseritropoyesis, 
  disgranulopoyesis, has_disgranulopoyesis,
  micromegacariocitos from patient p join analytical_data ad on 
  p.id = ad.patient_id join 
  morphological_data md on ad.id = md.morphologycal_data_id where register_number
  in (%s) and analysis_time = 1045 order by register_number", 
  paste0(all_num_web, collapse = ",")))

displasias<-data.frame(lapply(displasias, as.numeric))
displasias[displasias == 560]<-"Sí"
displasias[displasias == 550]<-"No"
displasias[displasias == 555]<-NA

displasias_1 <- displasias %>% mutate(my_diseritro = diseritropoyesis >= 10) %>% 
  mutate(my_disgranulopoyesis = disgranulopoyesis >= 10) %>% 
  mutate(my_megacariocitos = micromegacariocitos >= 10) %>% 
  rowwise() %>% mutate(suma_displasias = sum(my_megacariocitos, my_disgranulopoyesis,
                                             my_diseritro, na.rm = T)) %>% 
  mutate(multilinea = suma_displasias>1)

displasias_1$register_number<-as.numeric(displasias_1$register_number)

final_displasias<-merge(df_num_web, displasias_1, by.x = "num_web",
                        by.y = "register_number", all.x = T, all.Y = F)

#multilineage dysplasia? chequear displasias y ver si hay más de dos
#CREO que solucionado
#bone marrow fibrosis SOLUCIONADO

#### citogenetical_data: ISCN, number of metaphases ####

cariotype<-dbGetQuery(con,
                      sprintf("select register_number, cariotype_description, 
    analysis_date, ipss_cytogenetic, ipssr_cytogenetic from patient p join analytical_data ad
    on p.id = ad.patient_id join citogenetical_data cd on 
    cd.cito_data_id = ad.id where analysis_time = 1045 and
    register_number in (%s) ORDER BY register_number", 
                              paste0(all_num_web, collapse = ",")))

cariotype$register_number<-as.numeric(cariotype$register_number)

cariotype2<-dbGetQuery(con,
                       sprintf("select distinct on (register_number) register_number, cariotype_description, 
    analysis_date, ipss_cytogenetic, ipssr_cytogenetic from patient p join analytical_data ad
    on p.id = ad.patient_id join citogenetical_data cd on 
    cd.cito_data_id = ad.id where
    register_number in (%s) ORDER BY register_number", 
                               paste0(setdiff(all_num_web,cariotype$register_number), collapse = ",")))

cariotype2$register_number<-as.numeric(cariotype2$register_number)

cariotype <- rbind(cariotype,cariotype2)

cariotype$ipss_cytogenetic[cariotype$ipss_cytogenetic == 36119] <- 0
cariotype$ipss_cytogenetic[cariotype$ipss_cytogenetic == 36120] <- 0.5
cariotype$ipss_cytogenetic[cariotype$ipss_cytogenetic == 36121] <- 1
cariotype$ipssr_cytogenetic[cariotype$ipssr_cytogenetic == 36122] <- 0
cariotype$ipssr_cytogenetic[cariotype$ipssr_cytogenetic == 36123] <- 1
cariotype$ipssr_cytogenetic[cariotype$ipssr_cytogenetic == 36124] <- 2
cariotype$ipssr_cytogenetic[cariotype$ipssr_cytogenetic == 36125] <- 3
cariotype$ipssr_cytogenetic[cariotype$ipssr_cytogenetic == 36126] <- 4

cariotype$ipss_cytogenetic[cariotype$cariotype_description=="46,XY[20]" & is.na(cariotype$ipss_cytogenetic)] <- 0
cariotype$ipss_cytogenetic[cariotype$cariotype_description=="46,XX[20]" & is.na(cariotype$ipss_cytogenetic)] <- 0
cariotype$ipssr_cytogenetic[cariotype$cariotype_description=="46,XY[20]" & is.na(cariotype$ipssr_cytogenetic)] <- 1
cariotype$ipssr_cytogenetic[cariotype$cariotype_description=="46,XX[20]" & is.na(cariotype$ipssr_cytogenetic)] <- 1

cariotype$ipss_cytogenetic[grepl("(.*,.*){3,}", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipss_cytogenetic)] <- 1
cariotype$ipss_cytogenetic[grepl("-7|del\\(7q\\)", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipss_cytogenetic)] <- 1
cariotype$ipss_cytogenetic[grepl("-Y|del\\(5q\\)|del\\(20q\\)", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipss_cytogenetic)] <- 0
cariotype$ipss_cytogenetic[is.na(cariotype$ipss_cytogenetic)] <- 0


cariotype$ipssr_cytogenetic[grepl("(.*,.*){4,}", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 4
cariotype$ipssr_cytogenetic[grepl("(.*,.*){3}", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 3
cariotype$ipssr_cytogenetic[grepl("-7|inv\\(3\\)|t\\(3q\\)|del\\(3q\\)", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 3
cariotype$ipssr_cytogenetic[grepl("\\+8|del\\(7q\\)|i\\(17q\\)|\\+19", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 2
cariotype$ipssr_cytogenetic[grepl("del\\(20q\\)|del\\(5q\\)|del\\(12p\\)", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 1
cariotype$ipssr_cytogenetic[grepl("-Y|del\\(11q\\)", cariotype$cariotype_description) == TRUE & is.na(cariotype$ipssr_cytogenetic)] <- 0
cariotype$ipssr_cytogenetic[is.na(cariotype$ipssr_cytogenetic)] <- 1



final_cariotype<-merge(df_num_web, cariotype, by.x = "num_web",
                       by.y = "register_number", all.x = T, all.Y = F)

#### comorbidity ####

comorbidity <- dbGetQuery(con, sprintf(
  "select register_number, cd.* from patient p
  join comorbidity_data cd on p.id = cd.id where register_number in (%s) ORDER BY register_number", 
  paste0(all_num_web, collapse = ",")))

presente<-function(row){
  return(sum(as.numeric(row) == 18174, na.rm = T))
}

c2 <- comorbidity[,c(1,7:21)]
x <- sapply(c2, presente)

#### current status ####

estado_actual<-dbGetQuery(con, 
                          sprintf("with preseleccion as (select * from patient where register_number 
                      in (%s)) select register_number, assessment_date, death_reason,
                      status, disease_state_death from preseleccion pr join current_status cs
                      on pr.id = cs.id order by register_number",
                                  paste0(all_num_web, 
                                         collapse = ","))) 

estado_actual$register_number<-as.numeric(estado_actual$register_number)

estado_actual<-estado_actual %>% 
  merge(data.frame(num_web = all_num_web), 
        all.x = F, all.y = T, 
        by.x = "register_number", 
        by.y = "num_web")

codigos_status<-dbGetQuery(con, "select id, code from list_item where
                           id in (190, 195, 200)")

estado_actual2<-estado_actual %>% merge(codigos_status, by.x = "status", 
                                        by.y = "id", all.x = T, all.y = F) %>% 
  arrange(register_number)

estado_actual2$code<-gsub("ESTADO_PACIENTE_","", estado_actual2$code)

#### progresion ####

progresion<-dbGetQuery(con,
                       sprintf("select distinct(register_number), who2017_type_id, code as 
  tipo_progresion, progression, progression_date, evol_date from evolution_data ev
  join patient p on ev.patient_id = p.id join list_item li on ev.who2017_type_id
  = li.id where 
  register_number IN (%s) and who2017_type_id = 1034 order by 
  register_number;",
                               paste0(all_num_web, collapse = ","))) 

progresion$register_number<-as.numeric(progresion$register_number)

progresion<-progresion %>% merge(data.frame(num_web = all_num_web), 
                                 all.x = F, all.y = T, 
                                 by.x = "register_number", 
                                 by.y = "num_web") %>% 
  group_by(register_number) %>% arrange(desc(evol_date), .by_group = T) %>% 
  filter(row_number() == 1)

#### treatment ####

#iron chelation therapy es tratamiento quelante

trasplante <- dbGetQuery(con, sprintf(
  "select register_number, transplanted, tph_date from patient 
  join tph_data on patient.id = tph_data.id where register_number
  in (%s) order by register_number", paste0(all_num_web, collapse = ","))) %>% 
  filter(transplanted == TRUE) %>% select(-transplanted) %>% 
  mutate(tratamiento = "TPH") %>% rename(start_date = tph_date) %>% 
  select(register_number, tratamiento, start_date)  %>% 
  mutate(treatment_line = NA) %>% 
  mutate(end_date = NA)

trasplante$register_number <- as.numeric(trasplante$register_number)

transfusion <- dbGetQuery(con, sprintf(
  "select register_number, transfusion_date from patient p join 
  transfusion t on p.id = t.patient_id where register_number in (%s)
  order by register_number", paste0(all_num_web, collapse = ","))) %>% 
  mutate(tratamiento = "Transfusion") %>% 
  rename(start_date = transfusion_date) %>% 
  select(register_number, tratamiento, start_date) %>% 
  mutate(treatment_line = NA) %>% 
  mutate(end_date = NA)

transfusion$register_number <- as.numeric(transfusion$register_number)

#test <- dbGetQuery(con, "select * from treatment_data")
tratamiento<-dbGetQuery(con, sprintf(
  "select register_number, drug_name_id, support_drug_name_id,
   start_date, end_date, treatment_line, description
   from patient p join treatment_data td on p.id = td.patient_id
   left join treatment_drug_data tdd on td.id = tdd.treatment_data_id where
   register_number in (%s) order by register_number", 
  paste0(all_num_web, collapse = ","))) 


tratamiento$register_number<-as.numeric(tratamiento$register_number)

tratamiento2<-dbGetQuery(con, 
                         sprintf("select register_number, drug_name_id, support_drug_name_id,
          td.start_date, end_date, treatment_line, description
          from patient p join treatment_data td on p.id = td.patient_id
          join treatment_cycle_data tcd on td.id = tcd.treatment_data_id left join
          treatment_drug_data tdd2 on tcd.id = tdd2.cycle_data_id where
          register_number in (%s) order by register_number", 
                                 paste0(all_num_web, collapse = ","))) 

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

#View(lista_tratamientos)

t5 <- merge(t4, lista_tratamientos, by.x = "drug_id", by.y = "id",
            all.x = T, all.Y = F)

t5$tratamiento <- gsub("NOM_FARMACO_", "", t5$tratamiento)
t5$tratamiento <- gsub("NOMBRE_FARMACO_", "", t5$tratamiento)
t5$tratamiento <- ifelse(!is.na(t5$description),paste0(t5$tratamiento, "--", t5$description),t5$tratamiento)

t6 <- t5[,-1] %>% filter(!is.na(tratamiento)) %>% group_by(register_number) %>% 
  arrange(start_date, treatment_line, .by_group = T)

t6 <- t6 %>% filter(!is.na(tratamiento)) %>% group_by(register_number,treatment_line) %>%
  mutate(tratamiento=paste(tratamiento,collapse = "+")) %>%
  filter(!duplicated(tratamiento)) %>% ungroup(treatment_line) %>%
  arrange(start_date, treatment_line, .by_group = T)


t6$start_date[t6$start_date == "2999-12-12"] <- NA
t6$end_date[t6$end_date == "2999-12-12"] <- NA

# if (nrow(trasplante) > 0) {
#   t6_2 <- rbind(t6, trasplante)
# }else{
#   t6_2 <- t6
# }
# 
# if (nrow(transfusion) > 0) {
# t6_3 <- rbind(t6_2, transfusion) %>% arrange(start_date, treatment_line,
#                                              .by_group = T)
# }else{
#   t6_3 <- t6_2
# }



#View(t6)
#he quitado t6_3 de la linea de abajo porque añade cosas que no nos interesan (transfusión y trasplante con lineas de tto)

t7 <- t6 %>% select(register_number, tratamiento, start_date) %>% 
  summarise_all(function(x) x[1:4])

t8 <- t7 %>% unique() %>% group_by(register_number) %>% arrange(start_date, 
                                                                .by_group = T) %>% 
  summarise_all(function(x) x[1:4]) %>% rename(top4 = tratamiento)

t9 <- t8 %>% mutate(n=paste0("L", 1:n())) %>% 
  pivot_wider(names_from = n, values_from = c(top4, start_date))



colnames(t9)[2:9]<-c("Tratamiento_1", "Tratamiento_2", "Tratamiento_3",
                     "Tratamiento_4", "Fecha_tratamiento_1", 
                     "Fecha_tratamiento_2", "Fecha_tratamiento_3",
                     "Fecha_tratamiento_4")

t9<-t9[,c(1,2,6,3,7,4,8,5,9)]

#t9[t9 %in% c()] <- "HMA"
#t9[t9 %in% c()] <- "ESA"

t_quelante <- t5 %>% filter(tratamiento == "TRATAMIENTO_QUELANTE") %>% 
  select(register_number, start_date) %>% mutate(quelante = T)

final_quelante <- merge(df_num_web, t_quelante, by.x = "num_web",
                        by.y = "register_number", all.x = T, all.Y = F) %>% 
  group_by(num_web) %>% arrange(start_date, .by_group = T) %>% 
  filter(row_number() == 1)

final_t <- merge(df_num_web, t9, by.x = "num_web",
                 by.y = "register_number", all.x = T, all.y = F)

#### Intento de añadir mutaciones ####

mut3$register_number <- as.numeric(mut3$register_number)

mut4 <- mut3[,-1] %>% filter(!has_mutations == 0) %>% group_by(register_number) %>% 
  arrange(register_number, .by_group = T)

mut5 <- mut4 %>% filter(!is.na(vaf)) %>% select(register_number, gen, vaf) %>% 
  group_by(register_number) %>% summarise_all(function(x) x[1:9]) 

mut6 <- mut5 %>% mutate(n=paste0("_", 1:n())) %>% 
  pivot_wider(names_from = n, values_from = c(gen, vaf)) %>% arrange(register_number, .by_group = T)

mut6<-mut6[,c(1,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9,18,10,19)]

final_mut <- merge(df_num_web, mut6, by.x="num_web" ,by.y="register_number", all.x = T, all.y = F)

#### FINAL ####


final_df <- data.frame(patient_id = person$register_number,
                       nip = nhc$nip,
                       iniciales = nhc$identification,
                       center = person$hospital,
                       country = person$country,
                       gender = person$gender,
                       birthdate = person$birthdate,
                       diagnosis_date = diagnosis$diagnosis_date,
                       who_diagnosis_date = diagnosis$diagnosis_date,
                       who2017 = diagnosis$who2017,
                       who2022 = diagnosis$who2022,
                       icc2022 = diagnosis$icc2022,
                       dd = '',
                       ipss = final_pronostico$ipss,
                       ipssr = final_pronostico$ipssr,
                       ipss_blast = as.numeric(final_morphological$BM_BLAST_RS),
                       ipssr_blast = as.numeric(final_morphological$BM_BLAST_S),
                       ipss_cito = as.numeric(final_cariotype$ipss_cytogenetic),
                       ipssr_cito = as.numeric(final_cariotype$ipssr_cytogenetic),
                       ipss_hemo = as.numeric(final_hemogram$R),
                       ipssr_hb = as.numeric(final_hemogram$HB_RS),
                       ipssr_pmn = as.numeric(final_hemogram$ANC_RS),
                       ipssr_plq = as.numeric(final_hemogram$PLT_RS),
                       bone_marrow_blasts = final_morphological$blastos_mo,
                       periferal_blast = final_hemogram$blasts_percent,
                       age_diagnosis = diagnosis$edad_diagnostico,
                       date_hematological_collection = final_hemogram$analysis_date,
                       therapy_related = NA,
                       leuko = final_hemogram$leukocytes,
                       neutro = round(as.numeric(final_hemogram$neutrophils),2),
                       monocytes = round(as.numeric(final_hemogram$monocytes),2),
                       hemoglobin = final_hemogram$hemoglobin,
                       platelets = final_hemogram$platelets,
                       ldh = final_biochemistry$ldh,
                       ferritin = final_biochemistry$ferritin,
                       ring_sideroblasts = final_morphological$sideroblastos_anillados,
                       multilineage_dysplasia = final_displasias$multilinea, 
                       fibrosis_bm = final_morphological$fibrosis,
                       citogenetics = final_cariotype$cariotype_description,
                       metaphases = NA,
                       fish = NA,
                       comorbidity = "No",
                       type_comorbidity = NA,
                       date1 = final_t$Fecha_tratamiento_1,
                       treat1 = final_t$Tratamiento_1,
                       date2 = final_t$Fecha_tratamiento_2,
                       treat2 = final_t$Tratamiento_2,
                       date3 = final_t$Fecha_tratamiento_3,
                       treat3 = final_t$Tratamiento_3,
                       date4 = final_t$Fecha_tratamiento_4,
                       treat4 = final_t$Tratamiento_4,
                       gen1 = final_mut$gen__1,
                       vaf1 = final_mut$vaf__1,
                       gen2 = final_mut$gen__2,
                       vaf2 = final_mut$vaf__2,
                       gen3 = final_mut$gen__3,
                       vaf3 = final_mut$vaf__3,
                       gen4 = final_mut$gen__4,
                       vaf4 = final_mut$vaf__4,
                       gen5 = final_mut$gen__5,
                       vaf5 = final_mut$vaf__5,
                       gen6 = final_mut$gen__6,
                       vaf6 = final_mut$vaf__6,
                       gen7 = final_mut$gen__7,
                       vaf7 = final_mut$vaf__7,
                       gen8 = final_mut$gen__8,
                       vaf8 = final_mut$vaf__8,
                       gen9 = final_mut$gen__9,
                       vaf9 = final_mut$vaf__9,
                       iron_chelation = final_quelante$quelante,
                       iron_chelation_date = final_quelante$start_date,
                       evolution = !is.na(progresion$who2017_type_id),
                       evolution_date = progresion$evol_date,
                       outcome = estado_actual2$code
)

final_df$who2017[final_df$who2017=="WHO2017_LMMCX" & final_df$bone_marrow_blasts >= 10 ] <- "WHO2017_LMMC_2"
final_df$who2017[final_df$who2017=="WHO2017_LMMCX" & final_df$periferal_blast >= 5 ] <- "WHO2017_LMMC_2"
final_df$who2017[final_df$who2017=="WHO2017_LMMCX" & final_df$bone_marrow_blasts >= 5 ] <- "WHO2017_LMMC_1"
final_df$who2017[final_df$who2017=="WHO2017_LMMCX" & final_df$periferal_blast >= 2 ] <- "WHO2017_LMMC_1"
final_df$who2017[final_df$who2017=="WHO2017_LMMCX"] <- "WHO2017_LMMC_0"
final_df$who2022[final_df$who2022=="WHO2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_0"] <- "WHO2022_LMMC_1"
final_df$who2022[final_df$who2022=="WHO2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_1"] <- "WHO2022_LMMC_1"
final_df$who2022[final_df$who2022=="WHO2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_2"] <- "WHO2022_LMMC_2"
final_df$icc2022[final_df$icc2022=="ICC2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_0"] <- "ICC2022_LMMC_1"
final_df$icc2022[final_df$icc2022=="ICC2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_1"] <- "ICC2022_LMMC_1"
final_df$icc2022[final_df$icc2022=="ICC2022_LMMCX" & final_df$who2017 == "WHO2017_LMMC_2"] <- "ICC2022_LMMC_2"
final_df$dd[grepl("LMA",final_df$who2022)==TRUE] <- 1
final_df$dd[grepl("LMA",final_df$who2022)==FALSE] <- 2

final_df$ipss[is.na(final_df$ipss) & !is.na(final_df$ipss_blast) & !is.na(final_df$ipss_hemo) & !is.na(final_df$ipss_cito)] <- final_df$ipss_blast + final_df$ipss_cito + final_df$ipss_hemo
final_df$ipssr[is.na(final_df$ipssr) & !is.na(final_df$ipssr_blast) & !is.na(final_df$ipssr_cito) & !is.na(final_df$ipssr_hb) & !is.na(final_df$ipssr_pmn) & !is.na(final_df$ipssr_plq)] <- final_df$ipssr_blast + final_df$ipssr_cito + final_df$ipssr_hb + final_df$ipssr_pmn + final_df$ipssr_plq


comorbidity %>% filter_all(any_vars(. == 18174)) %>% .$register_number %>% 
  as.numeric() -> z
comorbidity %>% filter_all(any_vars(. == 18174)) %>% write.table("z.csv", 
                                                                 sep = ";",
                                                                 row.names = F,
                                                                 fileEncoding = "utf8")

#comorbidity %>% filter_all(any_vars(. == 18174))

final_df$comorbidity[final_df$patient_id %in% 
                       z] <- "Yes"

final_df2 <- as.data.frame(lapply(final_df, as.character))
final_df2$multilineage_dysplasia[is.na(final_df2$multilineage_dysplasia)] <- "No"
final_df2$fibrosis_bm[is.na(final_df2$fibrosis_bm)] <- "No"
final_df2[final_df2 == "FALSE"] <- "No"
final_df2[final_df2 == "TRUE"] <- "Yes"
final_df2[is.na(final_df2)] <- ""

final_df2$outcome[final_df2$outcome == "MUERTO"] <- "Dead"
final_df2$outcome[final_df2$outcome == "VIVO"] <- "Alive"
final_df2$outcome[final_df2$outcome == "PERDIDO"] <- "Lost"

bad <- read.delim(file.path("errores.txt"), header=FALSE)
bad <- as.vector(as.matrix(bad))
final_df2 <- preclean(final_df2,grep("citogenetics",colnames(final_df2)),bad)

final_df2$ipss[final_df2$ipss == ""] <- 999
final_df2$ipssr[final_df2$ipssr == ""] <- 999


#write_xlsx(final_df2,"ngs.xlsx")