#' ---------------------------
#'
#' Purpose of script:
#'
#'  Process clinical data used for training IPSS molecular (v2 after talking with Irene and Teresa)
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.0 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
clinical_raw <- read_tsv("./data/IPSSMol/df_clinical.tsv")

#' ---------------------------

# Edit variables
## Cytogenetics
#' Remove number of cell used for cytogenetics estimation
cytogenetics <- gsub("\\[.*\\]", "", clinical_raw$CYTOGENETICS)

#' Expand cytogenetic variables. Create a column per alteration and mark 
#' whether each individual had this alteration
cyto_list <- strsplit(cytogenetics, ",")

#' Simplify events
#' Select only frequent events reported by MDS foundation (https://www.mds-foundation.org/ipss-r-calculator/)
#' - LOY (-y)
#' - del(11q)
#' - del(5q)
#' - del(12p)
#' - del(20q)
#' - del(7q)
#' - +8
#' - +19
#' - i(17q)
#' - mono7
#' - inv(3)/t(3q)/del(3q)

simplify_deletion_addition <- function(annotation) {
  # Extract chromosome number and arm from the annotation
  extracted <- gsub("\\)\\(", "", annotation)
  extracted <- gsub("p.*$", "p)", extracted)
  extracted <- gsub("q.*$", "q)", extracted)

  # Return the simplified annotation
  return(extracted)
}

simplify_inversion <- function(annotation) {
  # Extract chromosome number and arm from the annotation
  out <-  regmatches(annotation, regexpr("inv\\(\\d+\\)", annotation))

  # Return the simplified annotation
  return(out)
}

simplify_translocations <- function(annotation) {
  # Extract chromosome number and arm from the annotation
  out <- gsub("t\\((.*);(.*)\\)\\(([pq]).*;([pq]).*\\)", "t(\\1\\3);t(\\2\\4)", annotation)


  # Return the simplified annotation
  return(out)
}



cyto_list_mod <- lapply(cyto_list, function(vec) {

    ## Remove underscores at the beggining
    vec <- gsub("^_", "", vec)
    dels <- vec[grep("^del", vec)]
    dup <- vec[grep("^dup", vec)]
    trans <- vec[grep("^t", vec)]
    inv <- vec[grep("^inv", vec)]
    iso <- vec[grep("^i\\(", vec)]
    other <- vec[!vec %in% c(dels, iso,  dup, trans, inv)]

    trans_out <- simplify_translocations(trans)
    c(other, simplify_deletion_addition(c(dels, dup, iso)), 
     unlist(strsplit(trans_out, ";")), simplify_inversion(inv) )
})
cyto_values <- unique(unlist(cyto_list_mod))

#' Create table with events
add_vars <- function(vec, events){

  out_vec <- rep(0, length(events))
  if(length(vec) == 0 && is.na(vec)){
    return(out_vec)
  }
  out_vec[events %in% vec] <- 1
  out_vec
}

events <- c("-y", "del(11q)", "del(5q)", "del(12p)", "del(20q)", 
  "del(7q)", "+8", "+19", "-7", "i(17q)", "inv(3)", "del(3q)", "t(3q)", "del(17p)")

cyto_tab <- sapply(cyto_list_mod, add_vars, events = events) %>% 
  t() %>% 
  as_tibble() 
colnames(cyto_tab) <- events

cyto_tab <- mutate(cyto_tab, "chr3ab" = `inv(3)`  + `del(3q)` + `t(3q)`) %>%
  select(-c(`inv(3)`, `del(3q)`, `t(3q)`))


#' Group less frequent events in the same category
#' Include a column to indicate if the patient had any aberration
#' If karyotype is complex, return maximum number of aberration (17)
createExtraVars <- function(vec) {
    if (length(vec) == 1 && is.na(vec)){
      return(c(0, 0))
    }
    filt_vec <- vec[!vec %in% c("46", "xy", "xx", "y", "x")]
    any_aberration <- 1
    if (length(filt_vec) == 0){
        any_aberration <- 0
    } 
    if (length(filt_vec) == 1 && filt_vec == "complex"){
      return(c(1, 17))
    }
    c(any_aberration, length(filt_vec))
}
cyto_tab_extra <- sapply(cyto_list_mod, createExtraVars) %>% 
  t() %>% 
  as_tibble() 
colnames(cyto_tab_extra) <- c("any_aberration", "N_aberrations")

cyto_tab_full <- cbind(cyto_tab, cyto_tab_extra) %>%
    as_tibble() %>%
    mutate(ID = clinical_raw$ID)

#' Adjust column names
colnames(cyto_tab_full) <- gsub("-", "mono", colnames(cyto_tab_full) )
colnames(cyto_tab_full) <- gsub("\\+", "triso", colnames(cyto_tab_full) )
colnames(cyto_tab_full) <- gsub("\\(|\\)", "", colnames(cyto_tab_full) )


## Copy Number
### Loss
loss_list <- strsplit(clinical_raw$CNACS_chrarm_loss, ",")
loss_events <- unique(unlist(loss_list))
loss_events <- loss_events[!is.na(loss_events)] ## Remove NAs

#' Select values with > 10 events
loss_tab <- table(unlist(loss_list))
loss_freq <- names(loss_tab[loss_tab > 10])

minor_val <- function(vec, events){

    if (length(vec) == 1 && is.na(vec)){
        return(0)
    } else if (all(vec %in% events)){
        return(0)
    } else {
        return(1)
    }
}

loss_mat <- sapply(loss_list, add_vars, events = loss_freq) %>% 
  t() %>% 
  as_tibble() 
colnames(loss_mat) <- paste0("loss_", loss_freq)

loss_mat <- loss_mat %>%
  mutate(loss_rare = sapply(loss_list, minor_val, events = loss_freq)) %>%
  mutate(loss_any = ifelse(rowSums(.) >= 1, 1, 0),
    ID = clinical_raw$ID)


### Gains
gain_list <- strsplit(clinical_raw$CNACS_chrarm_gain, ",")
gain_events <- unique(unlist(gain_list))
gain_events <- gain_events[!is.na(gain_events)] ## Remove NAs

#' Select values with > 10 events
gain_tab <- table(unlist(gain_list))
gain_freq <- names(gain_tab[gain_tab > 10])

gain_mat <- sapply(gain_list, add_vars, events = gain_freq) %>% 
  t() %>% 
  as_tibble() 
colnames(gain_mat) <- paste0("gain_", gain_freq)

gain_mat <- gain_mat %>%
  mutate(gain_rare = sapply(gain_list, minor_val, events = gain_freq)) %>%
  mutate(gain_any = ifelse(rowSums(.) >= 1, 1, 0),
    ID = clinical_raw$ID)


### UPD
upd_list <- strsplit(clinical_raw$CNACS_chrarm_upd, ",")
upd_events <- unique(unlist(upd_list))
upd_events <- upd_events[!is.na(upd_events)] ## Remove NAs

#' Select values with > 10 events
upd_tab <- table(unlist(upd_list))
upd_freq <- names(upd_tab[upd_tab > 10])

upd_mat <- sapply(upd_list, add_vars, events = upd_freq) %>% 
  t() %>% 
  as_tibble() 
colnames(upd_mat) <- paste0("upd_", upd_freq)

upd_mat <- upd_mat %>%
  mutate(upd_rare = sapply(upd_list, minor_val, events = upd_freq)) %>%
  mutate(upd_any = ifelse(rowSums(.) >= 1, 1, 0),
    ID = clinical_raw$ID)

# Data transform
#' Create levels for different variables
#' Round clinical variables 
#' Consider that samples with missing RINGED_SIDEROBLASTS have 0 sideroblasts.
clinical <- mutate(clinical_raw, 
    CYTO_IPSSR = factor(CYTO_IPSSR, levels = c("Very-Poor", "Poor", "Int", "Good", "Very-Good")),
    IPSSR = factor(IPSSR, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSRA = factor(IPSSRA, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
    OS_YEARS = ifelse(OS_YEARS == 0, 0.0001, OS_YEARS), ## Remove 0s in survival
    BM_BLAST = round(BM_BLAST, 0),
    PB_BLAST = round(PB_BLAST, 0),
    WBC = round(WBC, 0),
    ANC = round(ANC, 0), 
    MONOCYTES = round(MONOCYTES, 0),
    HB = round(HB, 0), 
    PLT = round(PLT, 0),
    RINGED_SIDEROBLASTS = ifelse(is.na(RINGED_SIDEROBLASTS), 0, RINGED_SIDEROBLASTS)
)

#' Merge all data.frames
clinical_all <- left_join(clinical, cyto_tab_full, by = "ID") %>%
    left_join(loss_mat, by = "ID") %>%
    left_join(gain_mat, by = "ID") %>%
    left_join(upd_mat, by = "ID")
save(clinical_all, file = "results/preprocess/clinical_preproc.Rdata")

