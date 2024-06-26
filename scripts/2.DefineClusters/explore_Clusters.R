#' ---------------------------
#'
#' Purpose of script:
#'
#' Explore clusters obtained from clinical variables
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.1 R
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
load("results/clustering/clinical_clustering_6.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(clin_clust6, levels = 1:6))

#' ---------------------------
# Explore clustering features 
#' ---------------------------

## Sex
#' ---------------------------
table(clinical_clust$clust, clinical_clust$SEX)
hma_prop <- prop.table(table(clinical_clust$clust, clinical_all$SEX), margin = 1)
hma_prop_df <- tibble(cluster = rownames(hma_prop), propF = hma_prop[, 1])
chisq.test(table(clinical_clust$clust, clinical_all$SEX))
hma_test <- lapply(levels(clinical_clust$clust), function(x) {
    test <- fisher.test(table(clinical_clust$clust == x, clinical_all$SEX))
    c(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
} ) %>%
    Reduce(rbind, .) %>%
    as_tibble() %>%
    mutate(cluster = as.character(cluster))

png("figures/clinclust_Sex.png")
left_join(hma_prop_df, hma_test, by = "cluster") %>%
    mutate(Significant = ifelse(pval < 0.05/14, "Significant", "Non-significant"),
        cluster = factor(cluster, levels = 1:14)) %>%
    ggplot(aes(y = propF*100, x = cluster, fill = Significant)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        geom_hline(yintercept = mean(clinical_all$SEX == "F")*100,
            linetype = "dashed") +
            ylab("Proportion of Females")
dev.off()
## Mayor mujeres: 5
## Mayor hombres: 2

## Age
#' ---------------------------
tapply(clinical_clust$AGE, clinical_clust$clust, summary)

png("figures/clinclust_age.png")
ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_all$AGE, na.rm = TRUE),
        linetype = "dashed")
dev.off()
summary(lm(AGE ~ clust, clinical_clust))

# Los pacientes del cluster 1 son más jóvenes

## IPSSM
#' ---------------------------
table(clinical_clust$clust, clinical_clust$IPSSM)
prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1)
chisq.test(table(clinical_clust$clust, clinical_clust$IPSSM))

png("figures/clinclust_ipssm.png")
table(clinical_clust$clust, clinical_clust$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/clinclust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

#' Los clusters tienen niveles de riesgo muy diferentes. 
#' Los pacientes de alto riesgo están en los clusteres 1 y 2.



## IPSSM Score
#' ---------------------------
tapply(clinical_clust$IPSSM_SCORE, clinical_clust$clust, summary)

png("figures/clinclust_ipssm_score.png")
ggplot(clinical_clust, aes(y = IPSSM_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

#' Las diferencias que veíamos usando la variable categórica se refuerzan
#' al considerar la variable continua. Para algunos clusteres (e.g.2 y 5), 
#' la variabilidad de los scores es bastante pequeña, sugieriendo que tienen
#' un riesgo similar.

## IPSSM Score vs IPSSRA
#' ---------------------------
png("figures/clinclust_ipssm_ipssra_score.png")
ggplot(clinical_clust, aes(x = IPSSRA_SCORE, y = IPSSM_SCORE)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~ clust)
dev.off()

sapply(levels(clinical_clust$clust), function(cl){
    a <- subset(clinical_clust, clust == cl)
    cor(a$IPSSRA_SCORE, a$IPSSM_SCORE, use = "complete")
})

## WHO Classification
#' ---------------------------
table(clinical_clust$clust, clinical_clust$WHO_2016)

png("figures/clinclust_who.png")
clust_who_tab <- table(clinical_clust$clust, clinical_clust$WHO_2016) %>%
    as.matrix()
clust_who_tab <- clust_who_tab[, !colnames(clust_who_tab) %in% c("MDS-RS-SLD/MLD", "aCML", "other")]
pheatmap(clust_who_tab, display_numbers = TRUE)
dev.off()

png("figures/clinclust_who_prop.png")
pheatmap(prop.table(clust_who_tab, margin = 2)*100, 
    display_numbers = TRUE)
dev.off()

png("figures/clinclust_who_prop2.png")
pheatmap(prop.table(clust_who_tab, margin = 1)*100, 
    display_numbers = TRUE)
dev.off()

#' Se ve una relación entre lus clusteres y los subptiopos de la OMS.
#' MDS RS: cluster 5
#' Del 5q: cluster 5
#' MDS MPNA: cluster 5
#' MDS EB2: cl2
#' CMML: cl 4 y 6
#' MDS MLD: cl 3.
#' 
#' Cl2: mayoría EB2


## Kariotipos
#' ---------------------------
kar_events <- c("monoy", "del11q", "del5q", "del12p",
                "del20q", "del7q", "triso8", "triso19", "mono7", 
                "i17q", "chr3ab", "complex")
clinical_clust$complex <- factor(clinical_all$complex, levels = c("non-complex", "complex"))

names(kar_events) <- kar_events
kar_tabs <- lapply(kar_events, function(x) table(clinical_clust[[x]], clinical_clust$clust))
kar_test <- lapply(kar_events, function(event) {
    lapply(levels(clinical_clust$clust), function(x) {
        test <- fisher.test(table(clinical_clust$clust == x, clinical_clust[[event]]))
        list(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
    } )     %>%
    Reduce(rbind, .) %>%
    as_tibble()  %>%
    mutate(cluster = factor(cluster, levels = 1:14), 
        OR = as.numeric(OR),
        pval = as.numeric(pval))
})
lapply(kar_test,function(x) arrange(x, OR))

kar_test_df <- Reduce(rbind, kar_test) %>%
    mutate(event = rep(names(kar_test), sapply(kar_test, nrow)),
        ORmod = ifelse(OR == 0, 0.18, OR),
        logOR = log(ORmod))

kar_test_m <- kar_test_df %>%
    select(cluster, event, logOR) %>%
    pivot_wider(names_from = event, values_from = logOR)
kar_test_mat <- as.matrix(kar_test_m[, -1])
rownames(kar_test_mat) <- kar_test_m$cluster

png("figures/clinclust_karevents_OR.png")
pheatmap(kar_test_mat, display_numbers = TRUE)
dev.off()

png("figures/clinclust_karevents_prop.png")
kar_tabs_prop <- sapply(kar_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = kar_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/clinclust_karevents_prop2.png")
kar_tabs_prop2 <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = kar_tabs_prop2, display_numbers = TRUE)
dev.off()

#' Se ven algunas correlaciones entre los clusters y los eventos kariotipicos

## Variables clinicas
#' ---------------------------
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clinical_clust, aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw()
})
png("figures/clinclust_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()

# Se observan bastantes diferencias a nivel de variables clinicas entre las muestras.



#' ---------------------------
# Train a decision tree 
#' ---------------------------
clinical_tree <- select(clinical_clust, c(clin_vars, "clust")) %>%
    mutate(logPLT = log(PLT),
            rat_BM_HB = BM_BLAST/HB,
            rat_WBC_ANC = WBC/ANC,
            rat_WBC_HB = WBC/HB,
            rat_WBC_PLT = WBC/logPLT) 

tree_mod <- rpart(clust ~ . , data = clinical_tree, method = "class", 
 control = rpart.control(cp = 0.005))

predictions <- predict(tree_mod, clinical_tree, type = "class")
clinical_tree$pred <- predictions

table(clinical_tree$pred, clinical_tree$clust)
sum(diag(table(predictions, clinical_tree$clust)))/length(predictions)

clin_tree_complete <- clinical_tree[complete.cases(clinical_tree), ]
table(clin_tree_complete$pred, clin_tree_complete$clust)
mean(clin_tree_complete$pred == clin_tree_complete$clust)

png("figures/classification_tree.png", width = 1400, height = 800)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()
#' ---------------------------
# Effect of clusters on survival 
#' ---------------------------

## Comparison with survival
#' ---------------------------
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_clust)

clinical_clust$clust <- relevel(clinical_clust$clust, "5") ## Lowest scores
png("figures/clinclust_survival.png")
os_clust <- survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()


## Se ven las diferencias que observábamos con los IPSSM scores


## Interacciones
#' ---------------------------
#' 
#' ## Sex
#' ---------------------------
cox_hma_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, clinical_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*hma + AGE + SEX, clinical_clust)

clinical_clust %>%
    subset(clust %in% c("4", "9", "11", "12", "7", "8", "14")) %>%
    mutate(clust = ifelse(clust %in% c("4", "9", "11"), "base", "other")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*clust + AGE, data = .)


hma_int_coef <- cox_hma_int$coef
hma_int_coef <- hma_int_coef[grep(":SEXM", names(hma_int_coef))]
hma_int_coef <- c("clust5:SEXM" = 0, hma_int_coef)
hma_int_coef <- gsub("clust", "", names(sort(hma_int_coef)))
hma_int_coef <- gsub(":SEXM", "", hma_int_coef)

png("figures/Surv_clustSex2.png", height = 800, width = 1200)
a <- lapply(hma_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    mutate(superClust = ifelse(clust %in% c(1, 4), "equal", "male")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*superClust + AGE, data = .)


#' ## HMA
#' ---------------------------

## Only high risk patients
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*IPSSM + SEX + AGE, data = clinical_clust)
coxph(Surv(OS_YEARS,OS_STATUS) ~ lenalidomid*IPSSM + SEX + AGE, data = clinical_clust)

filter(clinical_clust, IPSSM %in% c("High", "Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    group_by(hma, clust, IPSSM) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = hma, values_from = n) %>%
    mutate(n = `0` + `1`) %>%
    arrange(IPSSM, desc(n)) %>%
    print(n = Inf)

cox_hma_int_h <- clinical_clust %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

hma_int_coef_h <- cox_hma_int_h$coef
hma_int_coef_h <- hma_int_coef_h[grep("hma:", names(hma_int_coef_h))]
hma_int_coef_h <- c("hma:clust5" = 0, hma_int_coef_h)
hma_int_coef_h <- gsub("hma:clust", "", names(sort(hma_int_coef_h)))
png("figures/Surv_clusthma_high.png", width = 800)
a <- lapply(hma_int_coef_h, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

# clinical_clust %>%
#     subset(!clust %in% c(4, 9, 11, 12, 14, 10)) %>%
#     mutate(superClust = ifelse(clust %in% c(3, 5, 6), "low",
#         ifelse(clust %in% c(8, 13), "medium", "high"))) %>%
#     coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
clinical_clust %>%
    mutate(superClust = ifelse(clust %in% c(5, 4), "low",
        ifelse(clust %in% c(1, 3), "medium", "high"))) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(2, 6), "high", "low"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(2), "high", "low"))  %>%
    subset(clust != 6) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

names(hma_int_coef_h) <- hma_int_coef_h
lapply(hma_int_coef_h, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + SEX, 
        clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})
lapply(hma_int_coef_h, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + SEX + lenalidomid + chemotherapy +  transplant, 
        clinical_clust, subset = clust == x & IPSSM == "High" ) 
})

png("figures/Surv_clusthma_high_clustGroup.png")
cot <- clinical_clust %>%
    subset(IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(2, 6), "high", "low")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

## Very high
cox_hma_int_vh <- clinical_clust %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
hma_int_coef_vh <- cox_hma_int_vh$coef
hma_int_coef_vh <- hma_int_coef_vh[grep("hma:", names(hma_int_coef_vh))]
hma_int_coef_vh <- c("hma:clust5" = 0, hma_int_coef_vh)
hma_int_coef_vh <- gsub("hma:clust", "", names(sort(hma_int_coef_vh)))
png("figures/Surv_clusthma_veryhigh.png", width = 800)
a <- lapply(hma_int_coef_vh, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()


clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(2, 6), "high", "low"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
      mutate(superClust = ifelse(clust %in% c(2), "high", "low"))  %>%
    subset(clust != 6) %>%
     coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


names(hma_int_coef_vh) <- hma_int_coef_vh
lapply(hma_int_coef_vh, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + SEX, 
        clinical_clust, subset = clust == x & IPSSM == "Very-High"  & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})
lapply(hma_int_coef_vh, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + SEX + lenalidomid + chemotherapy +  transplant, 
        clinical_clust, subset = clust == x & IPSSM == "Very-High" ) 
})

png("figures/Surv_clusthma_veryhigh_clustGroup.png")
cot <- clinical_clust %>%
     subset(IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(2, 6), "high", "low")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

#' Parece que el cluster 2 es más efectivo a hma