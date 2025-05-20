#' ---------------------------
#'
#' Purpose of script:
#'
#' Explore meta-clusters obtained with all variables with umap
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
load("results/clustering/umap_meta_cluster.Rdata")
load("results/preprocess/ipssm_clinical_only.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(umap_pred)) %>%
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

png("figures/UMAPs/cluster_exploration/umapclust_Sex.png")
left_join(sex_prop_df, sex_test, by = "cluster") %>%
    mutate(Significant = ifelse(pval < 0.05/17, "Significant", "Non-significant")) %>%
        ggplot(aes(y = propF*100, x = cluster, fill = Significant)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        geom_hline(yintercept = mean(clinical_all$SEX == "F")*100,
            linetype = "dashed") +
            ylab("Proportion of Females")
dev.off()
#' Muchísimas diferencias en las proporciones entre hombres y mujeres. Hay 3 clusters donde 
#' solo ha ubicado mujeres (13, 56 y 57) y otros 4 predominantes de hombres (14, 16 y 4)


## Age
#' ---------------------------
tapply(clinical_clust$AGE, clinical_clust$clust, summary)

png("figures/UMAPs/cluster_exploration/umapclust_age.png")
ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_all$AGE, na.rm = TRUE),
        linetype = "dashed")
dev.off()
summary(lm(AGE ~ clust, clinical_clust))

#' Los pacientes del cluster 33 son los pacientes con MDS precoz. Hay algunas diferencias entre el resto, pero 
#' son menores.


## IPSSM
#' ---------------------------
table(clinical_clust$clust, clinical_clust$IPSSM)
prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1)
chisq.test(table(clinical_clust$clust, clinical_clust$IPSSM))

png("figures/UMAPs/cluster_exploration/umapclust_ipssm.png")
table(clinical_clust$clust, clinical_clust$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_ipssm_prop.png")
clist_ipssm_mat_prop <- prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1) %>%
    as.matrix()*100 
pheatmap(clist_ipssm_mat_prop, display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

#' Los clusters tienen niveles de riesgo muy diferentes. 
#' El cluster 6 contiene solo  pacientes de muy alto riesgo
#' Clusteres 1 y 5 tiene pacientes de alto y muy alto riesgo
#' Clusteres 14 y 33 tiene la mayor proporcion de pacientes de riesgo bajo o muy bajo
#' Otros clusteres también tienen más proporciones de pacientes de bajo riesgo.

## IPSSM Score
#' ---------------------------
tapply(clinical_clust$IPSSM_SCORE, clinical_clust$clust, summary)

png("figures/UMAPs/cluster_exploration/umapclust_ipssm_score.png")
ggplot(clinical_clust, aes(y = IPSSM_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

#' Las diferencias que veíamos usando la variable categórica se refuerzan
#' al considerar la variable continua. Para algunos clusteres (e.g.6, 14 o 33), 
#' la variabilidad de los scores es bastante pequeña, sugieriendo que tienen
#' un riesgo similar.


## IPSSRA Score
#' ---------------------------
tapply(clinical_clust$IPSSRA_SCORE, clinical_clust$clust, summary)

png("figures/UMAPs/cluster_exploration/umapclust_ipssra_score.png")
ggplot(clinical_clust, aes(y = IPSSRA_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

## IPSSM Score vs IPSSRA
#' ---------------------------
png("figures/UMAPs/cluster_exploration/umapclust_ipssm_ipssra_score.png")
ggplot(clinical_clust, aes(x = IPSSRA_SCORE, y = IPSSM_SCORE)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~ clust)
dev.off()

sapply(levels(clinical_clust$clust), function(cl){
    a <- subset(clinical_clust, clust == cl)
    cor(a$IPSSRA_SCORE, a$IPSSM_SCORE, use = "complete")
})
#' La correlación entre el IPSSM y el IPSSRA es muy variable dependiendo del cluster.
#' Tenemos clusters con una correlación muy alta (e.g. 4) que sugieren que las variables
#' que estaban incluidos en el IPSSRA (clinicas) son suficientemente predictivas del riesgo.
#' Por contra, otros clusters (eg. 57) tienen una correlación más pobre, sugiriendo que 
#' son necesarias las variables genéticas para afinar el riesgo.


## IPSSM Score vs IPSSM score clin only
#' ---------------------------
png("figures/UMAPs/cluster_exploration/umapclust_ipssm_ipssmclin_score.png")
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
#' Hay clusters con muy buena correlación (e.g. 98) y otros con una correlación muy pobre (e.g. 5). 


ipssm_tab <- table(molecular = clinical_clust$IPSSM, clinical = clinical_clust$IPSSMcat, clinical_clust$clust)

mask <- matrix(FALSE, 6, 6)
mask[3:6, 1] <- TRUE
mask[4:6, 2] <- TRUE
mask[5:6, 3] <- TRUE
mask[6, 4] <- TRUE

apply(ipssm_tab, 3, function(x) (sum(x) - sum(diag(x)))/sum(x)) %>% sort()
apply(ipssm_tab, 3, function(x) sum(x[mask])/sum(x)) %>% sort()

#' Hay algunos clusteres donde el IPSSM con los datos clínicos es parecido al 
#' IPSSM con datos moleculares también, como los clusters 57, 33 o 14. Estos pacientes
#' son de bajo riesgo y las diferencias entre pacientes en general se captan con
#' los datos clínicos.
#' Por otro lado, tenemos el cluster 6, donde todos los pacientes son de muy alto riesgo.
#' En este caso, ser del cluster 6 es más indicativo del riesgo que el IPSSM con datos clínicos.
#' 
#' Finalmente, tenemos otros casos donde el IPSSM con datos clínicos es muy poco predictivo. Por el 
#' ejemplo, en el cluster 1, 4 o 9, los datos moleculares incrementan el riesgo de muchos
#' pacientes.



coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*IPSSM_SCORE + AGE + SEX, clinical_clust)
## Revisar



## Raw survival
#' ---------------------------
png("figures/UMAPs/cluster_exploration/umapclust_raw_survival.png")
clinical_clust %>%
    filter(OS_STATUS == 1 | (OS_STATUS == 0 & OS_YEARS > 5)) %>%
    mutate(surv_years = ifelse(OS_YEARS > 5, 5, OS_YEARS)) %>%
    ggplot(aes(x = clust, y = surv_years)) +
        geom_boxplot() +
        theme_bw()
dev.off()


## WHO Classification
#' ---------------------------
table(clinical_clust$clust, clinical_clust$WHO_2016)

png("figures/UMAPs/cluster_exploration/umapclust_who.png")
clust_who_tab <- table(clinical_clust$clust, clinical_clust$WHO_2016) %>%
    as.matrix()
clust_who_tab <- clust_who_tab[, !colnames(clust_who_tab) %in% c("MDS-RS-SLD/MLD", "aCML", "other")]
pheatmap(clust_who_tab, display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_who_prop.png")
pheatmap(prop.table(clust_who_tab, margin = 2)*100, 
    display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_who_prop2.png")
pheatmap(prop.table(clust_who_tab, margin = 1)*100, 
    display_numbers = TRUE)
dev.off()

#' Se ve una relación entre lus clusteres y los subptipos de la OMS.
#' CMML: clusters 98 y 13
#' Del 5q: cluster 13
#' MDS/MPN-RS-T: clusters 18 y 57
#' MDS-EB2: Cluster 1 y 5

#' Cl1 y cl5: Mayoría MDS-EB2
#' Cl98: Mayoría CMML
#' Cl13: mayoría del5q
#' Cl33, 9, 14 y 56: MDS-MLD y MDS EB1


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
    mutate(OR = as.numeric(OR),
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

png("figures/UMAPs/cluster_exploration/umapclust_karevents_OR.png")
pheatmap(kar_test_mat, display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_karevents_prop.png")
kar_tabs_prop <- sapply(kar_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = kar_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_karevents_prop2.png")
kar_tabs_prop2 <- sapply(kar_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = kar_tabs_prop2, display_numbers = TRUE)
dev.off()

#' Se ven algunas correlaciones entre los clusters y los eventos kariotipicos
#' El cluster 6 está enriquecido en kariotipos complejos, mientras que
#' el cluster 13 está enriquecido en del5q.
#' El cluster 6 tiene una gran proporción de muestras con eventos kariotípicos, 
#' desde kariotipos complejos a triso19, del12p, del7q o mono7.
#' Mono 7 también tiene una proporción relevante en cl6.
#' Del5q se concentra en cl13, mientas que monoy en cl 4, 16 y 18.

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
    mutate(OR = as.numeric(OR),
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

png("figures/UMAPs/cluster_exploration/umapclust_mutations_OR.png", width = 1000, height = 800)
pheatmap(mut_test_mat, display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_mutation_prop.png", width = 1000, height = 800)
mut_tabs_prop <- sapply(mut_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = mut_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_mutation_prop2.png", width = 1000, height = 800)
mut_tabs_prop2 <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = mut_tabs_prop2, display_numbers = TRUE)
dev.off()

#' Mutacions de TP53 se concentran en el cluster 6. 
#' Mutaciones de STAG2, ML_PTD, DDX41 y IDH1 se concentran en el cluster 1. 
#' Mutaciones de STAG2, ML_PTD, ZSRS2, KRAS, U2AF1, ETV6, BCOR, IDH2, ASXL1 y RUNX1 se concentran en el cluster 4. 
#' Clusters 17 y 98 tiene una gran proporción de muestras con mutaciones en SRSF2 y TET2. 
#' Clusters 57, 16 y 18 tienen mayor proporción de SF3B1.


## Variables clinicas
#' ---------------------------
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clinical_clust, aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw()
})
png("figures/UMAPs/cluster_exploration/umapclust_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()

# Se observan bastantes diferencias a nivel de variables clinicas entre las muestras.

#' ---------------------------
# Effect of clusters on survival 
#' ---------------------------

## Comparison with survival
#' ---------------------------
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_clust)

clinical_clust$clust <- relevel(clinical_clust$clust, "33") ## Lowest scores
png("figures/UMAPs/cluster_exploration/umapclust_survival.png")
os_clust <- survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()

png("figures/UMAPs/cluster_exploration/umapclust_survival2.png", width = 1500)
a <- lapply(list(c(33, 14, 17, 57), c(13, 16, 18, 98), c(1, 4, 5, 6)), function(x){
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
    subset(!clust %in% c(14, 16, 4, 56, 57)) %>% ## Remove clusters with high sex-imbalance
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, data = .)

sex_int_coef <- cox_sex_int$coef
sex_int_coef <- sex_int_coef[grep(":SEXM", names(sex_int_coef))]
sex_int_coef <- c("clust33:SEXM" = 0, sex_int_coef)
sex_int_coef <- gsub("clust", "", names(sort(sex_int_coef)))
sex_int_coef <- gsub(":SEXM", "", sex_int_coef)

png("figures/UMAPs/cluster_exploration/umap_surv_clustSex.png", height = 800, width = 1200)
a <- lapply(sex_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    filter(!clust %in% c(14, 16, 4, 56, 57)) %>% ## Remove clusters with high sex-imbalance
    mutate(superClust = ifelse(clust %in% c(17, 5, 98, 1, 18 ), "equal", "male")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*superClust + AGE, data = .)



#' ## HMA
#' ---------------------------

## Only high risk patients
coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*IPSSM + SEX + AGE, data = clinical_clust)

filter(clinical_clust, IPSSM %in% c("High", "Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    group_by(hma, clust, IPSSM) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = hma, values_from = n) %>%
    mutate(n = `0` + `1`) %>%
    arrange(IPSSM, desc(n)) %>%
    print(n = Inf)

cox_hma_int_h <- clinical_clust %>%
    filter(clust %in% c(4, 1, 9, 56, 5, 98)) %>% ## Min 5 per group
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

hma_int_coef_h <- cox_hma_int_h$coef
hma_int_coef_h <- hma_int_coef_h[grep("hma:", names(hma_int_coef_h))]
hma_int_coef_h <- c("hma:clust1" = 0, hma_int_coef_h)
hma_int_coef_h <- gsub("hma:clust", "", names(sort(hma_int_coef_h)))
png("figures/UMAPs/cluster_exploration/umap_surv_clusthma_high.png", width = 800)
a <- lapply(hma_int_coef_h, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    filter(clust %in% c(4, 1, 9, 56, 5, 98)) %>% ## Min 5 per group
    mutate(superClust = ifelse(clust %in% c(4, 9), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

names(hma_int_coef_h) <- hma_int_coef_h
lapply(hma_int_coef_h, function(x){
    print(x)
    if (x == 56){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE, 
        clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
 
    } else{
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + SEX, 
        clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
    }
})

png("figures/UMAPs/cluster_exploration/Surv_clusthma_high_clustGroup.png")
cot <- clinical_clust %>%
    subset(IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(4, 9), "low", "high")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

## Very high
cox_hma_int_vh <- clinical_clust %>%
    filter(clust %in% c(6, 1, 4, 5, 9)) %>% ## Min 5 per group
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
hma_int_coef_vh <- cox_hma_int_vh$coef
hma_int_coef_vh <- hma_int_coef_vh[grep("hma:", names(hma_int_coef_vh))]
hma_int_coef_vh <- c("hma:clust1" = 0, hma_int_coef_vh)
hma_int_coef_vh <- gsub("hma:clust", "", names(sort(hma_int_coef_vh)))
png("figures/UMAPs/cluster_exploration/umap_surv_clusthma_veryhigh.png", width = 800)
a <- lapply(hma_int_coef_vh, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()


clinical_clust %>%
       mutate(superClust = ifelse(clust %in% c(6, 4), "low", "high"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)



names(hma_int_coef_vh) <- hma_int_coef_vh
lapply(hma_int_coef_vh, function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE , 
        clinical_clust, subset = clust == x & IPSSM == "Very-High"  & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
})


png("figures/UMAPs/cluster_exploration/umap_surv_clusthma_veryhigh_clustGroup.png")
cot <- clinical_clust %>%
     subset(IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(6, 4), "low", "high")) 
survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ sC + hma, cot) %>%
    ggsurvplot(data = cot) 
dev.off()

clinical_clust %>%
    filter(clust %in% c(1, 4, 5, 6, 9)) %>%
       mutate(superClust = ifelse(clust %in% c(4), "low", "high"))  %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM %in% c("High", "Very-High") & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

#' Parece que el cluster 4 es menos efectivo a hma



clinical_clust %>%
   subset(!WHO_2016 %in% c("MDS-RS-SLD", "MDS-RS-SLD/MLD", "MDS-SLD")) %>%
     coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*WHO_2016 + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

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

dd %>%
   subset(!who %in% c("other", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
   mutate(who = droplevels(who)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*who + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
