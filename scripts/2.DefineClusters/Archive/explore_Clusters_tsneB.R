#' ---------------------------
#'
#' Purpose of script:
#'
#' Explore clusters obtained with all variables with tsne
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
load("results/clustering/tsne10b_clusters.Rdata")
load("results/preprocess/ipssm_clinical_only.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(clust10_meta, levels = c(1:17, "5m"))) %>%
    left_join(ipssm_clinonly, by = "ID") %>%
    mutate(IPSSMcat = factor(IPSSMcat, levels = c("Very Low", "Low", "Moderate Low", "Moderate High", "High", "Very High")))

#' ---------------------------
# Explore clustering features 
#' ---------------------------

## Sex
#' ---------------------------
table(clinical_clust$clust, clinical_clust$SEX)
sex_prop <- prop.table(table(clinical_clust$clust, clinical_all$SEX), margin = 1)
sex_prop_df <- tibble(cluster = rownames(sex_prop), propF = sex_prop[, 1])
chisq.test(table(clinical_clust$clust, clinical_all$SEX))
sex_test <- lapply(levels(clinical_clust$clust), function(x) {
    test <- fisher.test(table(clinical_clust$clust == x, clinical_all$SEX))
    c(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
} ) %>%
    Reduce(rbind, .) %>%
    as_tibble() %>%
    mutate(cluster = as.character(cluster),
            pval = as.numeric(pval),
            OR = as.numeric(OR))

png("figures/tsne_10b/tsneclust_Sex.png")
left_join(sex_prop_df, sex_test, by = "cluster") %>%
    mutate(Significant = ifelse(pval < 0.05/18, "Significant", "Non-significant")) %>%
    ggplot(aes(y = propF*100, x = cluster, fill = Significant)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        geom_hline(yintercept = mean(clinical_all$SEX == "F")*100,
            linetype = "dashed") +
            ylab("Proportion of Females")
dev.off()
#' Muchísimas diferencias en las proporciones entre hombres y mujeres. Hay 4 clusters donde 
#' solo ha ubicado mujeres (9, 11, 14, 16) y otros 4 predominantes de hombres (2, 4, 23 y 13)


## Age
#' ---------------------------
tapply(clinical_clust$AGE, clinical_clust$clust, summary)

png("figures/tsne_10b/tsneclust_age.png")
ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_all$AGE, na.rm = TRUE),
        linetype = "dashed")
dev.off()
summary(lm(AGE ~ clust, clinical_clust))

#' Los pacientes del cluster 17 son los pacientes con MDS precoz. Hay algunas diferencias entre el resto, pero 
#' son menores.


## IPSSM
#' ---------------------------
table(clinical_clust$clust, clinical_clust$IPSSM)
prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1)
chisq.test(table(clinical_clust$clust, clinical_clust$IPSSM))

png("figures/tsne_10b/tsneclust_ipssm.png")
table(clinical_clust$clust, clinical_clust$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/tsne_10b/tsneclust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

#' Los clusters tienen niveles de riesgo muy diferentes. 
#' El cluster 6 contiene solo  pacientes de muy alto riesgo
#' Clusteres 3 y 5 tiene pacientes de alto y muy alto riesgo
#' Clusteres 7 y 12 tiene la mayor proporcion de pacientes de riesgo bajo o muy bajo
#' Otros clusteres también tienen más proporciones de pacientes de bajo riesgo.

## IPSSM Score
#' ---------------------------
tapply(clinical_clust$IPSSM_SCORE, clinical_clust$clust, summary)

png("figures/tsne_10b/tsneclust_ipssm_score.png")
ggplot(clinical_clust, aes(y = IPSSM_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

#' Las diferencias que veíamos usando la variable categórica se refuerzan
#' al considerar la variable continua. Para algunos clusteres (e.g.2, 6, 14 o 16), 
#' la variabilidad de los scores es bastante pequeña, sugieriendo que tienen
#' un riesgo similar.


## IPSSM Score vs IPSSM score clin only
#' ---------------------------
png("figures/tsne_10b/tsneclust_ipssm_ipssmclin_score.png")
ggplot(clinical_clust, aes(x = IPSSMscore, y = IPSSM_SCORE)) +
        geom_point() +
        theme_bw() +
        xlab("IPSSM clinical only") +
        ylab("IPSSM molecular") +
        facet_wrap(~ clust)
dev.off()

lapply(levels(clinical_clust$clust), function(cl){
    a <- subset(clinical_clust, clust == cl)
    summary(lm(IPSSMscore ~ IPSSM_SCORE, a))
})


sapply(levels(clinical_clust$clust), function(cl){
    a <- subset(clinical_clust, clust == cl)
    cor(a$IPSSMscore, a$IPSSM_SCORE, use = "complete")
})
#' Hay bastantes diferencias en la correlación entre el IPSSM con o sin datos moleculares dependiendo del cluster. 
#' Hay clusters con muy buena correlación (e.g. 1) y otros con una correlación muy pobre (e.g. 4). 


ipssm_tab <- table(molecular = clinical_clust$IPSSM, clinical = clinical_clust$IPSSMcat, clinical_clust$clust)

mask <- matrix(FALSE, 6, 6)
mask[3:6, 1] <- TRUE
mask[4:6, 2] <- TRUE
mask[5:6, 3] <- TRUE
mask[6, 4] <- TRUE

apply(ipssm_tab, 3, function(x) (sum(x) - sum(diag(x)))/sum(x)) %>% sort()
apply(ipssm_tab, 3, function(x) sum(x[mask])/sum(x)) %>% sort()

#' Hay algunos clusteres donde el IPSSM con los datos clínicos es parecido al 
#' IPSSM con datos moleculares también, como los clusters 14, 16, 12, 7 o 2. Estos pacientes
#' son de bajo riesgo y las diferencias entre pacientes en general se captan con
#' los datos clínicos.
#' Por otro lado, tenemos el cluster 6, donde todos los pacientes son de muy alto riesgo.
#' En este caso, ser del cluster 6 es más indicativo del riesgo que el IPSSM con datos clínicos.
#' 
#' Finalmente, tenemos otros casos donde el IPSSM con datos clínicos es muy poco predictivo. Por el 
#' ejemplo, en el cluster 4, 5 o 3, los datos moleculares incrementan el riesgo de muchos
#' pacientes.


coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*IPSSM_SCORE + AGE + SEX, clinical_clust)
#' Finalmente, vemos que el IPSSM_SCORE (con datos moleculares) es mucho menos predictivo de la supervivencia
#' en los clusters 6 y 8, sugiriendo que estos clusteres tienen pacientes con un riesgo más 
#' homogéneo y que no depende tanto de otras mutaciones. 

## WHO Classification
#' ---------------------------
table(clinical_clust$clust, clinical_clust$WHO_2016)

png("figures/tsne_10b/tsneclust_who.png")
clust_who_tab <- table(clinical_clust$clust, clinical_clust$WHO_2016) %>%
    as.matrix()
clust_who_tab <- clust_who_tab[, !colnames(clust_who_tab) %in% c("MDS-RS-SLD/MLD", "aCML", "other")]
pheatmap(clust_who_tab, display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_who_prop.png")
pheatmap(prop.table(clust_who_tab, margin = 2)*100, 
    display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_who_prop2.png")
pheatmap(prop.table(clust_who_tab, margin = 1)*100, 
    display_numbers = TRUE)
dev.off()

#' Se ve una relación entre lus clusteres y los subptiopos de la OMS.
#' CMML: cluster 10
#' Del 5q: cluster 9
#' MDS/MPN-RS-T: clusters 12 y 14
#' MDS-EB2: Cluster 3 y 5

#' Cl3, cl5 y cl6: Mayoría MDS-EB2
#' Cl10: Mayoría CMML
#' Cl9: mayoría del5q
#' Cl4, 11, 7 y 8: MDS-MLD y MDS EB1


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
    mutate(cluster = factor(cluster, levels = 1:17), 
        OR = as.numeric(OR),
        pval = as.numeric(pval))
})
lapply(kar_test,function(x) arrange(x, OR))

kar_test_df <- Reduce(rbind, kar_test) %>%
    mutate(event = rep(names(kar_test), sapply(kar_test, nrow)),
        ORmod = ifelse(OR == 0, 0.0395, OR),
        logOR = log(ORmod))

kar_test_m <- kar_test_df %>%
    select(cluster, event, logOR) %>%
    pivot_wider(names_from = event, values_from = logOR)
kar_test_mat <- as.matrix(kar_test_m[, -1])
rownames(kar_test_mat) <- kar_test_m$cluster

png("figures/tsne_10b/tsneclust_karevents_OR.png")
pheatmap(kar_test_mat, display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_karevents_prop.png")
kar_tabs_prop <- sapply(kar_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = kar_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_karevents_prop2.png")
kar_tabs_prop2 <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = kar_tabs_prop2, display_numbers = TRUE)
dev.off()

#' Se ven algunas correlaciones entre los clusters y los eventos kariotipicos
#' El cluster 6 está enriquecido en kariotipos complejos, mientras que
#' el cluster 9 está enriquecido en del5q.
#' El cluster 6 tiene una gran proporción de muestras con eventos kariotípicos, 
#' desde kariotipos complejos a triso19, del12p, del7q o mono7.
#' Mono 7 también tiene una proporción relevante en cl8.
#' Del5q se concentra en cl9, mientas que monoy en cl 12, 2 y 20.



## Mutaciones
#' ---------------------------
mutation <- read_tsv("./data/IPSSMol/df_mut.tsv")
sel_muts <- colMeans(mutation[, -1]) > 0.02 ## Select mutations present in at least 2% of the samples
sel_muts <- colnames(mutation[, -1])[sel_muts]

clin_mut <- left_join(clinical_clust, select(mutation, ID, sel_muts), by = "ID")

names(sel_muts) <- sel_muts
mut_tabs <- lapply(sel_muts, function(x) table(clin_mut[[x]], clin_mut$clust))
mut_test <- lapply(sel_muts, function(event) {
    lapply(levels(clin_mut$clust), function(x) {
        test <- fisher.test(table(clin_mut$clust == x, clin_mut[[event]]))
        list(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
    } )     %>%
    Reduce(rbind, .) %>%
    as_tibble()  %>%
    mutate(cluster = factor(cluster, levels = 1:17), 
        OR = as.numeric(OR),
        pval = as.numeric(pval))
})
lapply(mut_test,function(x) arrange(x, OR))

mut_test_df <- Reduce(rbind, mut_test) %>%
    mutate(event = rep(names(mut_test), sapply(mut_test, nrow)),
        ORmod = ifelse(OR == 0, 0.0395, OR),
        logOR = log(ORmod))

mut_test_m <- mut_test_df %>%
    select(cluster, event, logOR) %>%
    pivot_wider(names_from = event, values_from = logOR)
mut_test_mat <- as.matrix(mut_test_m[, -1])
rownames(mut_test_mat) <- mut_test_m$cluster

png("figures/tsne_10b/tsneclust_mutations_OR.png", width = 1000, height = 800)
pheatmap(mut_test_mat, display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_mutation_prop.png", width = 1000, height = 800)
mut_tabs_prop <- sapply(mut_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = mut_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/tsne_10b/tsneclust_mutation_prop2.png", width = 1000, height = 800)
mut_tabs_prop2 <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = mut_tabs_prop2, display_numbers = TRUE)
dev.off()

## Variables clinicas
#' ---------------------------
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clinical_clust, aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw()
})
png("figures/tsne_10b/tsneclust_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()

# Se observan bastantes diferencias a nivel de variables clinicas entre las muestras.



#' ---------------------------
# Train a decision tree 
#' ---------------------------
clinical_tree <- select(clinical_clust, c(clin_vars, "clust", "AGE", "SEX", kar_events)) %>%
    mutate(logPLT = log(PLT),
            rat_BM_HB = BM_BLAST/HB,
            rat_WBC_ANC = WBC/ANC,
            rat_WBC_HB = WBC/HB,
            rat_WBC_PLT = WBC/logPLT) 

tree_mod <- rpart(clust ~ . , data = clinical_tree, method = "class", 
 control = rpart.control(cp = 0.001))

predictions <- predict(tree_mod, clinical_tree, type = "class")
clinical_tree$pred <- predictions

table(clinical_tree$pred, clinical_tree$clust)
sum(diag(table(predictions, clinical_tree$clust)))/length(predictions)

png("figures/tsne_10b/tsne_classification_missing.png")
pheatmap(table(clinical_tree$pred, clinical_tree$clust), cluster_rows = FALSE, 
    cluster_cols = FALSE, display_numbers = TRUE, legend = FALSE)
dev.off()


clin_tree_complete <- clinical_tree[complete.cases(clinical_tree), ]
table(clin_tree_complete$pred, clin_tree_complete$clust)
mean(clin_tree_complete$pred == clin_tree_complete$clust)

png("figures/tsne_10b/tsne_classification_complete.png")
pheatmap(table(clin_tree_complete$pred, clin_tree_complete$clust), cluster_rows = FALSE, 
    cluster_cols = FALSE, display_numbers = TRUE, legend = FALSE)
dev.off()



png("figures/tsne_10b/tsne_classification_tree.png", width = 5000, height = 3500)
plot(tree_mod)
text(tree_mod, use.n = TRUE)
dev.off()
#' ---------------------------
# Effect of clusters on survival 
#' ---------------------------

## Comparison with survival
#' ---------------------------
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_clust)

clinical_clust$clust <- relevel(clinical_clust$clust, "7") ## Lowest scores
png("figures/tsne_10b/tsneclust_survival.png")
os_clust <- survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()

png("figures/tsne_10b/tsneclust_survival2.png", width = 1500)
a <- lapply(list(c(7, 14, 2, 12, 16, 10), c(1, 4, 9, 11, 13, 15, 17), c(3, 5, "5m", 6, 8)), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust, subset = clust %in% x) %>%
    ggsurvplot(data = clinical_clust)
    p$plot
})
plot_grid(plotlist = a, nrow = 1)
dev.off()

## Se ven las diferencias que observábamos con los IPSSM scores


## Interacciones
#' ---------------------------
#' 
#' ## Sex
#' ---------------------------
cox_sex_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, clinical_clust, subset = clust != 14)
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*hma + AGE + SEX, clinical_clust)

cox_sex_int <- clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, data = .)



sex_int_coef <- cox_sex_int$coef
sex_int_coef <- sex_int_coef[grep(":SEXM", names(sex_int_coef))]
sex_int_coef <- c("clust7:SEXM" = 0, sex_int_coef)
sex_int_coef <- gsub("clust", "", names(sort(sex_int_coef)))
sex_int_coef <- gsub(":SEXM", "", sex_int_coef)

png("figures/tsne_10b/tsne_surv_clustSex.png", height = 800, width = 1200)
a <- lapply(sex_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(superClust = ifelse(clust %in% c(10, 12, 3, "5m", 1), "equal", "male")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*superClust + AGE, data = .)

clinical_clust %>%
    subset(clust %in% c("3", "5", "5m")) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*clust + AGE, data = .)


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
    filter(clust %in% c(8, 5, "5m", 3, 4, 10, 11, 13, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

hma_int_coef_h <- cox_hma_int_h$coef
hma_int_coef_h <- hma_int_coef_h[grep("hma:", names(hma_int_coef_h))]
hma_int_coef_h <- c("hma:clust2" = 0, hma_int_coef_h)
hma_int_coef_h <- gsub("hma:clust", "", names(sort(hma_int_coef_h)))
png("figures/tsne_10b/tsne_surv_clusthma_high.png", width = 800)
a <- lapply(hma_int_coef_h, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    filter(clust %in% c(8, 5, 3, "5m", 4, 10, 11, 13, 2)) %>%
    mutate(superClust = ifelse(clust %in% c(2, 3, 10, 8), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

names(hma_int_coef_h) <- hma_int_coef_h
lapply(hma_int_coef_h, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE, 
        clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})


png("figures/tsne_10b/Surv_clusthma_high_clustGroup.png")
cot <- clinical_clust %>%
    subset(IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(2, "5m"), "high", "low")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

## Very high
cox_hma_int_vh <- clinical_clust %>%
    filter(clust %in% c(6, 5, "5m", 3, 8, 11, 13)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
hma_int_coef_vh <- cox_hma_int_vh$coef
hma_int_coef_vh <- hma_int_coef_vh[grep("hma:", names(hma_int_coef_vh))]
hma_int_coef_vh <- c("hma:clust3" = 0, hma_int_coef_vh)
hma_int_coef_vh <- gsub("hma:clust", "", names(sort(hma_int_coef_vh)))
png("figures/tsne_10b/tsne_surv_clusthma_veryhigh.png", width = 800)
a <- lapply(hma_int_coef_vh, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()


clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(6, 3, 13), "low", "high"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)



names(hma_int_coef_vh) <- hma_int_coef_vh
lapply(hma_int_coef_vh, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE , 
        clinical_clust, subset = clust == x & IPSSM == "Very-High"  & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})


png("figures/tsne_10b/tsne_surv_clusthma_veryhigh_clustGroup.png")
cot <- clinical_clust %>%
     subset(IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(6, 3, 13), "low", "high")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(3), "low", "high"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM %in% c("High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust == 3 & IPSSM %in% c("High", "Very-High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
clinical_clust %>%
    subset(clust == 3 & IPSSM %in% c("High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust == 3 & IPSSM %in% c("Very-High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
    subset(clust == 5 & IPSSM %in% c("High", "Very-High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
clinical_clust %>%
    subset(clust == 5 & IPSSM %in% c("High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust == 5 & IPSSM %in% c("Very-High")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust %in% c(5, 3) & IPSSM %in% c("High", "Very-High")) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust %in% c(5, 3) & IPSSM %in% c("High")) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
    subset(clust %in% c(5, 3) & IPSSM %in% c("Very-High")) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
    subset(clust %in% c(5, 11, 3, 6) & IPSSM %in% c("High", "Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6) & IPSSM %in% c("High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6) & IPSSM %in% c("Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
    subset(clust %in% c(5, 11, 3, 6, 8, "5m") & IPSSM %in% c("High", "Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6, "5m"), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6, "5m") & IPSSM %in% c("High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6, "5m"), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6, 8, "5m") & IPSSM %in% c("Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6, "5m"), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

#' Parece que el cluster 3 es menos efectivo a hma
clinical_clust %>%
   subset(!WHO_2016 %in% c("MDS-RS-SLD", "MDS-RS-SLD/MLD", "MDS-SLD")) %>%
     coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*WHO_2016 + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)



## Comparison with AML transformation
#' ---------------------------
coxph(Surv(AMLt_YEARS, AMLt_STATUS) ~ clust + AGE + SEX, clinical_clust)


png("figures/tsne_10b/tsneclust_amlt.png")
os_clust <- survfit(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()

png("figures/tsne_10b/tsneclust_amlt2.png", width = 1500)
a <- lapply(list(c(16, 14, 12, 2), c(7, 10, 13, 17, 15, 4, 9, 11), c(1, 3, 5, "5m", 6, 8)), function(x){
p <-  survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ clust, clinical_clust, subset = clust %in% x) %>%
    ggsurvplot(data = clinical_clust)
    p$plot
})
plot_grid(plotlist = a, nrow = 1)
dev.off()

## Se ven las diferencias que observábamos con los IPSSM scores

#' ## Sex
#' ---------------------------
cox_sex_int_aml <- clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ clust*SEX + AGE, data = .)



sex_int_coef_aml <- cox_sex_int_aml$coef
sex_int_coef_aml <- sex_int_coef_aml[grep(":SEXM", names(sex_int_coef_aml))]
sex_int_coef_aml <- c("clust7:SEXM" = 0, sex_int_coef_aml)
sex_int_coef_aml <- gsub("clust", "", names(sort(sex_int_coef_aml)))
sex_int_coef_aml <- gsub(":SEXM", "", sex_int_coef_aml)

png("figures/tsne_10b/tsne_amlt_clustSex.png", height = 800, width = 1200)
a <- lapply(sex_int_coef_aml, function(x){
p <-  survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(superClust = ifelse(clust %in% c(10, 15), "0female",
    ifelse(clust %in% c(16, 17, 12, 1), "2male", "1equal"))) %>%
    coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ SEX*superClust + AGE, data = .)



#' ## HMA
#' ---------------------------
coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ hma + SEX + AGE, data = clinical_clust, subset = IPSSM %in% c("Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

filter(clinical_clust, IPSSM %in% c("High", "Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    group_by(hma, clust, IPSSM, AMLt_STATUS) %>%
    summarize(n = n()) %>%
    filter(!is.na(AMLt_STATUS)) %>%
    pivot_wider(names_from = hma, values_from = n) %>%
    mutate(n = `0` + `1`) %>%
    arrange(IPSSM, clust) %>%
    print(n = Inf)

## Very high
cox_hma_int_aml_vh <- clinical_clust %>%
    filter(clust %in% c(5, 6, 8, "5m")) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
hma_int_coef_aml_vh <- cox_hma_int_aml_vh$coef
hma_int_coef_aml_vh <- hma_int_coef_aml_vh[grep("hma:", names(hma_int_coef_aml_vh))]
hma_int_coef_aml_vh <- c("hma:clust5" = 0, hma_int_coef_aml_vh)
hma_int_coef_aml_vh <- gsub("hma:clust", "", names(sort(hma_int_coef_aml_vh)))
png("figures/tsne_10b/tsne_amlt_clusthma_veryhigh.png", width = 800)
a <- lapply(hma_int_coef_aml_vh, function(x){
p <-  survfit(formula = Surv(AMLt_YEARS,AMLt_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()


clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(5, "5m"), "low", "high"))  %>%
    coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)



names(hma_int_coef_vh) <- hma_int_coef_vh
lapply(hma_int_coef_vh, function(x){
    coxph(Surv(AMLt_YEARS,AMLt_STATUS) ~ hma + AGE + SEX, 
        clinical_clust, subset = clust == x & IPSSM == "Very-High"  & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})


## Comparison with LFS
#' ---------------------------
coxph(Surv(LFS_YEARS, LFS_STATUS) ~ clust + AGE + SEX, clinical_clust)


png("figures/tsne_10b/tsneclust_lfs.png")
os_clust <- survfit(Surv(LFS_YEARS,LFS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()

png("figures/tsne_10b/tsneclust_lfs2.png", width = 1500)
a <- lapply(list(c(14, 12, 2, 7), c(10, 17, 1, 9, 13, 15, 4, 11), c(3, 5, "5m", 6, 8)), function(x){
p <-  survfit(formula = Surv(LFS_YEARS,LFS_STATUS) ~ clust, clinical_clust, subset = clust %in% x) %>%
    ggsurvplot(data = clinical_clust)
    p$plot
})
plot_grid(plotlist = a, nrow = 1)
dev.off()

## Se ven las diferencias que observábamos con los IPSSM scores

#' ## Sex
#' ---------------------------
cox_sex_int_lfs <- clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(LFS_YEARS,LFS_STATUS) ~ clust*SEX + AGE, data = .)



sex_int_coef_lfs <- cox_sex_int_lfs$coef
sex_int_coef_lfs <- sex_int_coef_lfs[grep(":SEXM", names(sex_int_coef_lfs))]
sex_int_coef_lfs <- c("clust7:SEXM" = 0, sex_int_coef_lfs)
sex_int_coef_lfs <- gsub("clust", "", names(sort(sex_int_coef_lfs)))
sex_int_coef_lfs <- gsub(":SEXM", "", sex_int_coef_lfs)

png("figures/tsne_10b/tsne_lfs_clustSex.png", height = 800, width = 1200)
a <- lapply(sex_int_coef_lfs, function(x){
p <-  survfit(formula = Surv(LFS_YEARS,LFS_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    subset(!clust %in% c("4",  "9", "13", "11", "14")) %>% ## Remove clusters with high sex-imbalance
    mutate(superClust = ifelse(clust %in% c(10, 15, "5m", 12), "equal", "male")) %>%
    coxph(Surv(LFS_YEARS,LFS_STATUS) ~ SEX*superClust + AGE, data = .)


#' ## HMA
#' ---------------------------
coxph(Surv(LFS_YEARS,LFS_STATUS) ~ hma + SEX + AGE, data = clinical_clust, subset = IPSSM %in% c("Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

filter(clinical_clust, IPSSM %in% c("High", "Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    group_by(hma, clust, IPSSM, LFS_STATUS) %>%
    summarize(n = n()) %>%
    filter(!is.na(LFS_STATUS)) %>%
    pivot_wider(names_from = hma, values_from = n) %>%
    mutate(n = `0` + `1`) %>%
    arrange(IPSSM, clust) %>%
    print(n = Inf)


dd <- clinical_clust
dd$who <- dd$WHO_2016
dd$who[dd$who %in% c("MDS-SLD","MDS-MLD")] <- "MDS-SLD/MLD"
dd$who[dd$who %in% c("MDS-RS-SLD","MDS-RS-MLD")] <- "MDS-RS-SLD/MLD"
dd$who[dd$who %in% c("MDS-EB1","MDS-EB2")] <- "MDS-EB1/2"
dd$who[dd$WHO_2016=="aCML"] <- "other"
dd$who <- factor(as.vector(dd$who),levels=c("MDS-del5q",
                                             "MDS-RS-SLD/MLD","MDS-SLD/MLD",
                                              "MDS-EB1/2","MDS-U","CMML","MDS/MPN-RS-T",
                                              "MDS/MPN-U","other"))
dd %>%
   subset(!who %in% c("other", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
   mutate(who = droplevels(who)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*who + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

                                              "MDS/MPN-U","other"))
dd %>%
   subset(!who %in% c("other", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
   mutate(who = droplevels(who)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*who + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
