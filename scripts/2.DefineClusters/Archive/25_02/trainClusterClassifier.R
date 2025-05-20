#' ---------------------------
#'
#' Purpose of script:
#'
#' Train a classifier to get clusters in other datasets
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.3 R
#'
#' ---------------------------


# Load libraries and data
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)
library(cowplot)
library(rpart)

load("results/preprocess/clinical_preproc.Rdata")
load("results/tsne_perp10.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(clust10_tree, levels = 1:17))

#' ---------------------------
# Train a decision tree 
# Assume missing PB and monocytes are 0
#' ---------------------------
kar_events <- c("N_aberrations", "del5q")
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "logPLT")

data_tree <- mutate(clinical_clust, 
    logPLT = log(PLT),
    MONOCYTES = ifelse(is.na(MONOCYTES), 0, MONOCYTES),
    PB_BLAST = ifelse(is.na(PB_BLAST), 0, PB_BLAST)) %>%
    select(c(clin_vars, "clust", "AGE", "SEX", kar_events)) 
data_tree_comp <- data_tree[complete.cases(data_tree),]

clust_freq <- data_tree_comp %>%
    group_by(clust) %>%
    summarize(w = n()/nrow(data_tree_comp)) %>%
    pull(w)
names(clust_freq) <- 1:17

weights <- 1/clust_freq[data_tree_comp$clust]

tree_mod <- rpart(clust ~ . , weights = weights, data = data_tree_comp, method = "class")

predictions <- predict(tree_mod, data_tree_comp, type = "class")
acc_tab <- table(prediction = predictions, original = data_tree_comp$clust)
sum(diag(acc_tab))/length(predictions)


diag(acc_tab)/colSums(acc_tab)
diag(acc_tab)/rowSums(acc_tab)

png("figures/tsne_clust_classifier/classification_tree_raw.png", width = 2000, height = 1000)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()


data_tree_comp2 <- subset(data_tree_comp, !clust %in% c(3, 5, 6, 10)) %>%
    mutate(clust = droplevels(clust))
tree_mod2 <- rpart(clust ~ . , weights = 1/clust_freq[data_tree_comp2$clust], 
    data = data_tree_comp2, method = "class")


predictions2 <- predict(tree_mod2, data_tree_comp2, type = "class")
acc_tab2 <- table(prediction = predictions2, original = data_tree_comp2$clust)
sum(diag(acc_tab2))/length(predictions2)


diag(acc_tab2)/colSums(acc_tab2)
diag(acc_tab2)/rowSums(acc_tab2)

png("figures/tsne_clust_classifier/classification_tree_raw2.png", width = 2000, height = 1000)
plot(tree_mod2)
text(tree_mod2, use.n = TRUE)
dev.off()



## Differentiate cluster 2 from 4, 7, 15
clust2_sub <- subset(data_tree, clust %in% c(2, 4, 7, 15))


comp_vars <- c("WBC", "ANC", "HB", "logPLT", "AGE", "SEX")

pdf("figures/tsne_clust_classifier/clust2_diffs.pdf")
lapply(seq_len(length(comp_vars) - 1), function(i){
    lapply(seq(from = i + 1, to = length(comp_vars)), function(j){
    ggplot(clust2_sub, aes_string(x = comp_vars[i], y =  comp_vars[j], col = "clust")) +
        geom_point(position = "jitter")
    })
})
dev.off()

## Differentiate cluster 17 from 4, 11 and 15
clust17_sub <- subset(data_tree, clust %in% c(4, 11, 15, 17))

pdf("figures/tsne_clust_classifier/clust17_diffs.pdf")
lapply(seq_len(length(comp_vars) - 1), function(i){
    lapply(seq(from = i + 1, to = length(comp_vars)), function(j){
    ggplot(clust17_sub, aes_string(x = comp_vars[i], y =  comp_vars[j], col = "clust")) +
        geom_point(position = "jitter")
    })
})
dev.off()

### Try "manual classification"
data_tree_mini <- mutate(clinical_clust, 
    logPLT = log(PLT),
    MONOCYTES = ifelse(is.na(MONOCYTES), 0, MONOCYTES),
    PB_BLAST = ifelse(is.na(PB_BLAST), 0, PB_BLAST),
    WBC_HB = WBC/HB) %>%
    filter(complex == "non-complex" & del5q == 0 & MONOCYTES < 1.5 & BM_BLAST < 8.5) %>%
    select(c(clin_vars, "clust", "AGE", "SEX", WBC_HB)) 
data_tree_mini <- data_tree_mini[complete.cases(data_tree_mini),] %>%
    filter(!clust %in% c(3, 5, 6, 9, 10)) %>%
    mutate(clust = droplevels(clust))



tree_mini <- rpart(clust ~ . , weights = 1/clust_freq[data_tree_mini$clust],
    data = data_tree_mini, method = "class")
pred_mini <- predict(tree_mini, data_tree_mini, type = "class")
acc_mini_tab <- table(prediction = pred_mini, original = data_tree_mini$clust)
sum(diag(acc_mini_tab))/length(pred_mini)




## Test SVM
library(e1071)
svm_mod <- svm(clust ~ . , data = data_tree_comp, kernel = "linear")
pred_svm <- predict(svm_mod, data_tree_comp, type = "class")
table(pred_svm, data_tree_comp$clust)
sum(diag(table(pred_svm, data_tree_comp$clust)))

svm_mod <- svm(clust ~ . , data = data_tree_comp)
pred_svm <- predict(svm_mod, data_tree_comp, type = "class")
table(pred_svm, data_tree_comp$clust)
sum(diag(table(pred_svm, data_tree_comp$clust)))

svm_mod <- svm(clust ~ . , data = data_tree_comp, kernel = "polynomial")
pred_svm <- predict(svm_mod, data_tree_comp, type = "class")
table(pred_svm, data_tree_comp$clust)
sum(diag(table(pred_svm, data_tree_comp$clust)))


clinical_tree$pred <- predictions

table(clinical_tree$pred, clinical_tree$clust)
sum(diag(table(predictions, clinical_tree$clust)))/length(predictions)

png("figures/tsne_classification_missing.png")
pheatmap(table(clinical_tree$pred, clinical_tree$clust), cluster_rows = FALSE, 
    cluster_cols = FALSE, display_numbers = TRUE, legend = FALSE)
dev.off()


clin_tree_complete <- clinical_tree[complete.cases(clinical_tree), ]
table(clin_tree_complete$pred, clin_tree_complete$clust)
mean(clin_tree_complete$pred == clin_tree_complete$clust)

png("figures/tsne_classification_complete.png")
pheatmap(table(clin_tree_complete$pred, clin_tree_complete$clust), cluster_rows = FALSE, 
    cluster_cols = FALSE, display_numbers = TRUE, legend = FALSE)
dev.off()



png("figures/tsne_classification_tree.png", width = 5000, height = 3500)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()