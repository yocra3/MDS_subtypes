#' ---------------------------
#'
#' Purpose of script:
#'
#' Compute manually IPSSM score to make a fair comparison
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run --rm -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.6 R
#'
#' ---------------------------

## Load data and libraries
library(tidyverse)
library(ipssm)
library(survival)

clinical_raw <- read_tsv("./data/IPSSMol/df_clinical.tsv")
cyto <- read_table("data/IPSSMol/df_cna.tsv")
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
maf <- read.table("./data/IPSSMol/maf.tsv", header=T, sep="\t", stringsAsFactors = F)

clin_train <- read_table("results/gnn/preprocess/patient_variables.tsv")

clinical_train <- left_join(clinical_raw,  select(clin_train, ID, train), by = "ID") %>%
    left_join(cyto, by = "ID") %>%
    left_join(mutation, by = "ID") 

## Modify variables to match IPSSM names
clinical_train$del7_7q <- as.numeric(clinical_train$del7q | clinical_train$del7)
clinical_train$del17_17p <- as.numeric(clinical_train$del17p | clinical_train$del17)
clinical_train$complex <- car::recode(clinical_train$complex, "'complex'='1'; 'non-complex'='0'")
clinical_train$CYTO_IPSSR <- car::recode(clinical_train$CYTO_IPSSR, "'Very-Good'='Very Good' ; 'Int'='Intermediate' ; 'Very-Poor'='Very Poor'")
clinical_train$TP53mut <- as.vector(sapply(clinical_train$ID, function(x) sum(maf$ID==x & maf$GENE=="TP53")))
clinical_train$TP53mut <- car::recode(clinical_train$TP53mut, " '2'='2 or more' ; '3'='2 or more' ")
clinical_train$TP53loh <- 0 
clinical_train$TP53loh[which(clinical_train$chr17 %in% c("cnloh","del"))] <- 1
clinical_train$TP53maxvaf <- NA
clinical_train$FLT3[which(clinical_train$FLT3_ITD==1)] <- 1

## Process data for model
clin.process <- IPSSMprocess(clinical_train) %>%
    mutate(nRes2 = nRes2mean)

model_vars <- names(formals(IPSSMmain)$betaValues)
model_vars <- model_vars[model_vars != ""]
clinical_train_vars <- clin.process %>%
    select(all_of(model_vars), "AGE", "SEX", "train", "LFS_YEARS", "LFS_STATUS", "IPSSM_SCORE", "ID")  %>%
    as_tibble() %>%
    filter(!is.na(LFS_YEARS) & !is.na(LFS_STATUS))

cox_df <- clinical_train_vars %>% 
    select(-train, -IPSSM_SCORE, -ID)

all_cox_model <- coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ., data = cox_df)


## Train only on training set
cox_df_train <- clinical_train_vars %>% 
    filter(is.na(train) | train) %>%
    select(-train, -IPSSM_SCORE, -ID)

train_cox_model <- coxph(Surv(LFS_YEARS, LFS_STATUS) ~ ., data = cox_df_train)

clinical_train_vars_pred <- clinical_train_vars %>%
    mutate(AGE = 0, SEX = "F") %>%
    mutate(across(all_of(model_vars), function(x) x - mean(x, na.rm = TRUE))) %>%
    mutate(IPSSM_SCORE_new = log(predict(all_cox_model, newdata = ., type = "risk")),
           IPSSM_SCORE_train = log(predict(train_cox_model, newdata = ., type = "risk")))
cor(clinical_train_vars_pred$IPSSM_SCORE_new, clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs")
cor(clinical_train_vars_pred$IPSSM_SCORE_train, clinical_train_vars_pred$IPSSM_SCORE, use = "pairwise.complete.obs")

png("figures/gnn/recomputeIPSSM.png", width = 800, height = 600)
plot(clinical_train_vars_pred$IPSSM_SCORE_new, clinical_train_vars_pred$IPSSM_SCORE)
dev.off()


png("figures/gnn/recomputeIPSSM_train.png", width = 800, height = 600)
clinical_train_vars_pred %>%
    mutate(train = ifelse(is.na(train) | train, "Train", "Test")) %>%
    ggplot(aes(x = IPSSM_SCORE_train, y = IPSSM_SCORE)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(title = "Recomputed IPSSM (Train)",
         x = "Recomputed IPSSM (Train)",
         y = "Original IPSSM") +
    facet_grid(~ train) +
    theme_bw()
dev.off()

write.table(clinical_train_vars_pred, 
            file = "results/gnn/patient_recomputedIPSSM.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)