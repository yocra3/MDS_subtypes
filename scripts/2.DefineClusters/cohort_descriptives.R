#' ---------------------------
#'
#' Purpose of script:
#'
#'  Descriptives of the cohort
#' 
#' ---------------------------
#'
#' Notes:
#' Compute the descriptives of GESMD and IWS cohorts. Include patients with
#' classification and patients used for clustering.
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------


# Load libraries and data
library(tidyverse)
