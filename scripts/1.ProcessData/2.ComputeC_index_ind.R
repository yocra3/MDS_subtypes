#' ---------------------------
#'
#' Purpose of script:
#'
#'  Compute individual C-index for IPSSM 
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.5 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(tidyverse)
load("results/preprocess/clinical_preproc.Rdata")

score_mat <- clinical[, c("IPSSM_SCORE", "OS_YEARS", "OS_STATUS")]

samps <- nrow(score_mat)
rownames(score_mat) <- clinical$ID
risk_row <- matrix(score_mat$IPSSM_SCORE, nrow = samps, ncol = samps)
risk_col <- matrix(score_mat$IPSSM_SCORE, nrow = samps, ncol = samps, byrow = TRUE)
risk_bool <- risk_row > risk_col
risk_bool2 <- risk_row == risk_col

time_row <- matrix(score_mat$OS_YEARS, nrow = samps, ncol = samps)
time_col <- matrix(score_mat$OS_YEARS, nrow = samps, ncol = samps, byrow = TRUE)
time_bool <- time_row < time_col

event_row <- matrix(score_mat$OS_STATUS, nrow = samps, ncol = samps)
event_col <- matrix(score_mat$OS_STATUS, nrow = samps, ncol = samps, byrow = TRUE)
valid <- time_row < time_col & event_row == 1 | time_row > time_col & event_col == 1   

c_index_mat <- (risk_bool == time_bool)
c_index_mat[(risk_bool2 * valid) == 1] <- 0.5
c_index_mat[!valid] <- NA
rownames(c_index_mat) <- clinical$ID

ov_ind_c <- rowMeans(c_index_mat, na.rm = TRUE)

subtypes <- unique(clinical$consensus)
subtype_ind_c <- lapply(subtypes, function(subtype){

    sel_samps <- subset(clinical, consensus == subtype)$ID
    vec <- rownames(c_index_mat) %in% sel_samps
    mini_mat <- c_index_mat[vec, vec]
    ind_c <- rowMeans(mini_mat, na.rm = TRUE)
}) %>% unlist()

c_index_df <- data.frame(ID = names(ov_ind_c), OV_C_INDEX = ov_ind_c, SUB_C_INDEX = subtype_ind_c[names(ov_ind_c)])

clinical_c_index <- left_join(clinical, c_index_df, by = "ID")
save(clinical_c_index, file = "results/preprocess/clinical_c_index.Rdata")
