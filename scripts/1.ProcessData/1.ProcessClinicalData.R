#' ---------------------------
#'
#' Purpose of script:
#'
#'  Process clinical data used for training IPSS molecular.
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
cytogenetics <- gsub("\\[.*\\]", "", clinical$CYTOGENETICS)

#' Expand cytogenetic variables. Create a column per alteration and mark 
#' whether each individual had this alteration
cyto_list <- strsplit(cytogenetics, ",")

#' Simplify events
#' - Deletions (del), additions (add) and duplications (dup): Indicate just the 
#' chromosome and arm
#' - Derivatives (der): indicate the chromosome
#' - Translocations (t), dicentric (dic) and inversions (inv): only indicate
#'  there was a translocation, inversion or dicentric chromosome
#' - Marker chromosomes: group any marker chromosome indication
#' - Group rare monosomies (<44) and trisomies (>47).
#' 
simplify_deletion_addition <- function(annotation) {
  # Extract chromosome number and arm from the annotation
  extracted <- gsub("\\)\\(", "", annotation)
  extracted <- gsub("p.*$", "p)", extracted)
  extracted <- gsub("q.*$", "q)", extracted)

  # Return the simplified annotation
  return(extracted)
}

simplify_derivate <- function(annotation) {
  # Extract chromosome number and arm from the annotation
  out <-  regmatches(annotation, regexpr("der\\(\\d+\\)", annotation))

  # Return the simplified annotation
  return(out)
}
cyto_list_mod <- lapply(cyto_list, function(vec) {

    ## Remove underscores at the beggining
    vec <- gsub("^_", "", vec)
    dels <- vec[grep("^del", vec)]
    adds <- vec[grep("^add", vec)]
    dup <- vec[grep("^dup", vec)]
    der <- vec[grep("^der", vec)]
    trans <- vec[grep("^t", vec)]
    inv <- vec[grep("^inv", vec)]
    dic <- vec[grep("^dic", vec)]
    mar <- vec[grep("mar", vec)]
    other <- vec[!vec %in% c(dels, adds, der, dup, trans, inv, dic, mar)]

    ## Group monosomies and trisomies
    other[other %in% c("41-45", "41~44", "41~46", "42", "42-45", "42~45", 
        "43", "43-44", "43-45", "43~44", "43~46")] <- "<44"

    other[other %in% c("46~48", "47~48", "47~49", "47~55", "48",
     "4n", "178~182<7n>")] <- ">47"
    
    c(other, simplify_deletion_addition(c(dels,adds, dup)), 
    simplify_derivate(der), rep("trans", length(trans)), 
    rep("inv", length(inv)), rep("dic", length(dic)),
    rep("mar", length(mar))) 
})

cyto_values <- unique(unlist(cyto_list_mod))

#' Select values with > 10 events
cyto_vals_tab <- table(unlist(cyto_list_mod))
cyto_freq <- names(cyto_vals_tab[cyto_vals_tab > 10])

#' Remove non-descriptive values
#' Remove marker chromosomes and correct karyotypes
cyto_freq_sel <- cyto_freq[!cyto_freq %in% c("46", "xy", "xx", "y", "x")]

#' Create table with events
add_vars <- function(vec, events){

  out_vec <- rep(0, length(events))
  if(length(vec) == 0 && is.na(vec)){
    return(out_vec)
  }
  out_vec[events %in% vec] <- 1
  out_vec
}

cyto_tab <- sapply(cyto_list_mod, add_vars, events = cyto_freq_sel) %>% 
  t() %>% 
  as_tibble() 
colnames(cyto_tab) <- cyto_freq_sel

#' Group less frequent events in the same category
#' Include a column to indicate if the patient had any aberration
createExtraVars <- function(vec, events) {

    filt_vec <- vec[!vec %in% events]
    extra_del <- any(grepl("^del", filt_vec)) %>%
        as.numeric()
    extra_der <- any(grepl("^der", filt_vec)) %>%
        as.numeric()
    extra_add <- any(grepl("^add", filt_vec)) %>%
        as.numeric()
    extra_dup <- any(grepl("^dup", filt_vec)) %>%
        as.numeric()
    any_aberration <- 1
    if (length(vec) == 0 && is.na(vec)){
        any_aberration <- 0
    } 
    else if (all(vec %in% c("46", "xy", "xx", "y", "x"))){
        any_aberration <- 0
    }
    c(extra_del, extra_der, extra_add, extra_dup, any_aberration)
}
cyto_tab_extra <- sapply(cyto_list_mod, createExtraVars, events = cyto_freq_sel) %>% 
  t() %>% 
  as_tibble() 
colnames(cyto_tab_extra) <- c("del_rare", "der_rare", "add_rare", "dup_rare", "any_aberration")

cyto_tab_full <- cbind(cyto_tab, cyto_tab_extra) %>%
    as_tibble() %>%
    mutate(ID = clinical$ID,
           del_any = ifelse(rowSums(select(., starts_with("del"))) >= 1, 1, 0),
           der_any = ifelse(rowSums(select(., starts_with("der"))) >= 1, 1, 0),
           monosomies = ifelse(rowSums(select(., starts_with("-"), "44", "45", "<44")) >= 1, 1, 0),
           trisomies = ifelse(rowSums(select(., starts_with("+"), "47", ">47")) >= 1, 1, 0))

#' Adjust column names
colnames(cyto_tab_full) <- gsub("-", "mono", colnames(cyto_tab_full) )
colnames(cyto_tab_full) <- gsub("\\+", "triso", colnames(cyto_tab_full) )
colnames(cyto_tab_full) <- gsub("\\(|\\)", "", colnames(cyto_tab_full) )

colnames(cyto_tab_full)[grep("^[<>4]", colnames(cyto_tab_full))] <- 
    paste0("kar", colnames(cyto_tab_full)[grep("^[<>4]", colnames(cyto_tab_full))])

## Copy Number
### Loss
loss_list <- strsplit(clinical$CNACS_chrarm_loss, ",")
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
    ID = clinical$ID)


### Gains
gain_list <- strsplit(clinical$CNACS_chrarm_gain, ",")
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
    ID = clinical$ID)


### UPD
upd_list <- strsplit(clinical$CNACS_chrarm_upd, ",")
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
    ID = clinical$ID)

# Data transform
#' Create levels for different variables
clinical <- mutate(clinical_raw, 
    CYTO_IPSSR = factor(CYTO_IPSSR, levels = c("Very-Poor", "Poor", "Int", "Good", "Very-Good")),
    IPSSR = factor(IPSSR, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSRA = factor(IPSSRA, levels = c("Very-Low", "Low", "Int", "High", "Very-High")),
    IPSSM = factor(IPSSM, levels = c("Very-Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very-High")),
    OS_YEARS = ifelse(OS_YEARS == 0, 0.0001, OS_YEARS) ## Remove 0s in survival
)

#' Merge all data.frames
clinical <- left_join(clinical, cyto_tab_full, by = "ID") %>%
    left_join(loss_mat, by = "ID") %>%
    left_join(gain_mat, by = "ID") %>%
    left_join(upd_mat, by = "ID")
save(clinical_all, file = "results/preprocess/clinical_preproc.Rdata")
