library(tidyverse)

load("results/gesmd_data_all.Rdata")

## Filter patients
gesmd_valencia <- gesmd_data %>%
  mutate(ID = as.character(register_number),
         WHO_2016 = who2017) %>%
  filter(!is.na(TP53)) %>%
  filter(hospital %in% c("H. G. de Valencia", "H. C.  Valencia"))

tab <- gesmd_valencia %>%
  select(register_number, hospital, del5q) %>%
  mutate(register_number = as.character(register_number))

write.table(tab, file = "valencia_gesmd.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)