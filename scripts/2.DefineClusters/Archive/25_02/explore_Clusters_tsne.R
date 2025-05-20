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
load("results/tsne_perp10.Rdata")
load("results/preprocess/ipssm_clinical_only.Rdata")

clinical_clust <- mutate(clinical_all, 
    clust = factor(clust10_tree, levels = 1:17)) %>%
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

png("figures/tsneclust_Sex.png")
left_join(sex_prop_df, sex_test, by = "cluster") %>%
    mutate(Significant = ifelse(pval < 0.05/17, "Significant", "Non-significant"),
        cluster = factor(cluster, levels = 1:17)) %>%
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

png("figures/tsneclust_age.png")
ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_all$AGE, na.rm = TRUE),
        linetype = "dashed")
dev.off()
summary(lm(AGE ~ clust, clinical_clust))

age_plot <- ggplot(clinical_clust, aes(y = AGE, x = clust)) +
        geom_boxplot() +
        theme_bw() +
        geom_hline(yintercept = median(clinical_all$AGE, na.rm = TRUE),
        linetype = "dashed") +
        xlab("Cluster") +
        ylab("Age")




#' Los pacientes del cluster 17 son los pacientes con MDS precoz. Hay algunas diferencias entre el resto, pero 
#' son menores.


## IPSSM
#' ---------------------------
table(clinical_clust$clust, clinical_clust$IPSSM)
prop.table(table(clinical_clust$clust, clinical_clust$IPSSM), margin = 1)
chisq.test(table(clinical_clust$clust, clinical_clust$IPSSM))

png("figures/tsneclust_ipssm.png")
table(clinical_clust$clust, clinical_clust$IPSSM) %>%
    as.matrix() %>%
    pheatmap(display_numbers = TRUE, cluster_cols = FALSE)
dev.off()

png("figures/tsneclust_ipssm_prop.png")
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

png("figures/tsneclust_ipssm_score.png")
ggplot(clinical_clust, aes(y = IPSSM_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

#' Las diferencias que veíamos usando la variable categórica se refuerzan
#' al considerar la variable continua. Para algunos clusteres (e.g.2, 6, 14 o 16), 
#' la variabilidad de los scores es bastante pequeña, sugieriendo que tienen
#' un riesgo similar.


## IPSSRA Score
#' ---------------------------
tapply(clinical_clust$IPSSRA_SCORE, clinical_clust$clust, summary)

png("figures/tsneclust_ipssra_score.png")
ggplot(clinical_clust, aes(y = IPSSRA_SCORE, x = clust)) +
        geom_boxplot() +
        theme_bw()
dev.off()

## IPSSM Score vs IPSSRA
#' ---------------------------
png("figures/tsneclust_ipssm_ipssra_score.png")
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
#' Tenemos clusters con una correlación muy alta (e.g. 1 o 15) que sugieren que las variables
#' que estaban incluidos en el IPSSRA (clinicas) son suficientemente predictivas del riesgo.
#' Por contra, otros clusters (eg. 14 o 16) tienen una correlación más pobre, sugiriendo que 
#' son necesarias las variables genéticas para afinar el riesgo.


## IPSSM Score vs IPSSM score clin only
#' ---------------------------
png("figures/tsneclust_ipssm_ipssmclin_score.png")
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

png("figures/tsneclust_who.png")
clust_who_tab <- table(clinical_clust$clust, clinical_clust$WHO_2016) %>%
    as.matrix()
clust_who_tab <- clust_who_tab[, !colnames(clust_who_tab) %in% c("MDS-RS-SLD/MLD", "aCML", "other")]
pheatmap(clust_who_tab, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_who_prop.png")
pheatmap(prop.table(clust_who_tab, margin = 2)*100, 
    display_numbers = TRUE)
dev.off()

png("figures/tsneclust_who_prop2.png")
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

who_types <- colnames(clust_who_tab)
names(who_types) <- who_types
who_test <- lapply(who_types, function(who) {
    lapply(levels(clinical_clust$clust), function(x) {
        test <- fisher.test(table(clinical_clust$clust == x, clinical_clust$WHO_2016 == who))
        list(cluster = x, OR = as.numeric(test$estimate), pval = test$p.value)
    } )     %>%
    Reduce(rbind, .) %>%
    as_tibble()  %>%
    mutate(cluster = factor(cluster, levels = 1:17), 
        OR = as.numeric(OR),
        pval = as.numeric(pval))
})
lapply(who_test,function(x) arrange(x, OR))
lapply(who_test[c("MDS-SLD", "MDS-MLD")],function(x) arrange(x, OR))



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

png("figures/tsneclust_karevents_OR.png")
pheatmap(kar_test_mat, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_karevents_prop.png")
kar_tabs_prop <- sapply(kar_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = kar_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_karevents_prop2.png")
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

mutation_N <- mutate(mutation, 
    N_mutations = rowSums(mutation[, 2:127]),
    mutation0 = as.numeric(N_mutations == 0),
    mutation1 = as.numeric(N_mutations == 1),
    mutation2_3 = as.numeric(N_mutations %in% c(2, 3)),
    mutation_multi = as.numeric(N_mutations > 3))
    


clin_mut <- left_join(clinical_clust, select(mutation_N, ID, sel_muts, mutation0, mutation1, mutation2_3, mutation_multi), by = "ID")

sel_muts <- c(sel_muts, "mutation0", "mutation1", "mutation2_3", "mutation_multi")
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

png("figures/tsneclust_mutations_OR.png", width = 1000, height = 800)
pheatmap(mut_test_mat, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_mutation_prop.png", width = 1000, height = 800)
mut_tabs_prop <- sapply(mut_tabs, function(x) prop.table(x, margin = 1)[2, ])*100
pheatmap(mat = mut_tabs_prop, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_mutation_prop2.png", width = 1000, height = 800)
mut_tabs_prop2 <- sapply(mut_tabs, function(x) prop.table(x, margin = 2)[2, ])*100
pheatmap(mat = mut_tabs_prop2, display_numbers = TRUE)
dev.off()

## Compare with consensus classification
#' ---------------------------
clin_mut <- clin_mut %>%
    mutate(consensus = ifelse(TP53multi == 1 & BM_BLAST <= 20, "Mutated TP53",
        ifelse(del5q == 1 & del7q == 0 & BM_BLAST <= 5, "del5q",
            ifelse(SF3B1 > 0 & del7q == 0 & complex == "non-complex" & BM_BLAST <= 5, "mutated SF3B1",
                ifelse(BM_BLAST <= 5, "Low blasts",
                    ifelse(BM_BLAST > 10, "MDS-IB2",
                        ifelse(BM_BLAST > 5 & BM_BLAST <= 10, "MDS-IB1", "Other")))))))

table(clin_mut$clust, clin_mut$consensus)

png("figures/tsneclust_consensus.png")
clust_consensus_tab <- table(clin_mut$clust, clin_mut$consensus) %>%
    as.matrix()
pheatmap(clust_consensus_tab, display_numbers = TRUE)
dev.off()

png("figures/tsneclust_consensus_prop.png")
pheatmap(prop.table(clust_consensus_tab, margin = 2)*100, 
    display_numbers = TRUE)
dev.off()

png("figures/tsneclust_consensus_prop2.png")
pheatmap(prop.table(clust_consensus_tab, margin = 1)*100, 
    display_numbers = TRUE)
dev.off()

## Variables clinicas
#' ---------------------------
clin_vars <- c("BM_BLAST", "PB_BLAST", "WBC", "ANC", "MONOCYTES", "HB", "PLT")

clin_plots <- lapply(clin_vars, function(x){
   ggplot(clinical_clust, aes(y = .data[[x]], x = clust)) +
    geom_boxplot() +
    theme_bw()
})
png("figures/tsneclust_clinical_vars.png", width = 1000)
plot_grid(plotlist = clin_plots, nrow = 3)
dev.off()

# Se observan bastantes diferencias a nivel de variables clinicas entre las muestras.



#' ---------------------------
# Effect of clusters on survival 
#' ---------------------------

## Comparison with survival
#' ---------------------------
cox_surv <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust + AGE + SEX, clinical_clust)

clinical_clust$clust <- relevel(clinical_clust$clust, "7") ## Lowest scores
png("figures/tsneclust_survival.png")
os_clust <- survfit(Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust)
ggsurvplot(os_clust, data = clinical_clust) 
dev.off()


png("figures/tsneclust_survival2.png", width = 1500)
a <- lapply(list(c(7, 14, 2, 12, 16, 10), c(1, 4, 9, 11, 13, 15, 17), c(3, 5, 6, 8)), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ clust, clinical_clust, subset = clust %in% x) %>%
    ggsurvplot(data = clinical_clust)
    p$plot
})
plot_grid(plotlist = a, nrow = 1)
dev.off()

## Se ven las diferencias que observábamos con los IPSSM scores



## Compare C-index
clin_comp <- filter(clinical_clust, !is.na(OS_YEARS)) %>%
    mutate(clin_score = -IPSSMscore, 
    mol_score = -IPSSM_SCORE,
    ori_score = -IPSSRA_SCORE)

clust_model <- coxph(Surv(OS_YEARS, OS_STATUS) ~ clust, data = clin_comp)
clust_risk <- -predict(clust_model, type = "risk", data = clin_comp)

clin_who <- filter(clin_comp, !is.na(WHO_2016))
who_model <- coxph(Surv(OS_YEARS, OS_STATUS) ~ WHO_2016, data = clin_who)
who_risk <- -predict(who_model, type = "risk", data = clin_who)


concordance(Surv(OS_YEARS, OS_STATUS) ~ clust_risk, data = clin_comp)
# n= 2811 
# Concordance= 0.6527 se= 0.008182
# concordant discordant     tied.x     tied.y    tied.xy 
#    1269761     647127     122192        537         49 

concordance(Surv(OS_YEARS, OS_STATUS) ~ who_risk, data = clin_who)
# n= 2756 
# Concordance= 0.6406 se= 0.007968
# concordant discordant     tied.x     tied.y    tied.xy 
#    1123086     573515     258303        450        103 

clust_model2 <- coxph(Surv(OS_YEARS, OS_STATUS) ~ clust, data = clin_who)
clust_risk2 <- predict(clust_model2, type = "risk", data = clin_who)

con_clust <- survcomp::concordance.index(clust_risk2, clin_who$OS_YEARS, clin_who$OS_STATUS, method = "noether")
con_who <- survcomp::concordance.index(-who_risk, clin_who$OS_YEARS, clin_who$OS_STATUS, method="noether")
survcomp::cindex.comp(con_clust, con_who)

concordance(Surv(OS_YEARS, OS_STATUS) ~ ori_score, data = clin_comp)
# n=2555 (256 observations deleted due to missingness)
# Concordance= 0.7153 se= 0.007935
# concordant discordant     tied.x     tied.y    tied.xy 
#    1170454     465437       1727        454          0 

concordance(Surv(OS_YEARS, OS_STATUS) ~ clin_score, data = clin_comp)
# n=2609 (202 observations deleted due to missingness)
# Concordance= 0.7151 se= 0.007773
# concordant discordant     tied.x     tied.y    tied.xy 
#    1230357     488536       5449        480          0 
concordance(Surv(OS_YEARS, OS_STATUS) ~ mol_score, data = clin_comp)
# n=2596 (215 observations deleted due to missingness)
# Concordance= 0.7483 se= 0.007228
# concordant discordant     tied.x     tied.y    tied.xy 
#    1276943     428497       2941        478          0 


clin_full <- filter(clin_comp, !is.na(ori_score) & !is.na(clin_score) & !is.na(mol_score))
con_ipssr <- survcomp::concordance.index(-clin_full$ori_score, clin_full$OS_YEARS, clin_full$OS_STATUS, method = "noether")
con_clin <- survcomp::concordance.index(-clin_full$clin_score, clin_full$OS_YEARS, clin_full$OS_STATUS, method = "noether")
con_ipssm <- survcomp::concordance.index(-clin_full$mol_score, clin_full$OS_YEARS, clin_full$OS_STATUS, method = "noether")

survcomp::cindex.comp(con_ipssr, con_clin)
# $p.value
# [1] 0.5779894
survcomp::cindex.comp(con_ipssm, con_ipssr)
# $p.value
# [1] 8.329199e-14
survcomp::cindex.comp(con_ipssm, con_clin)
# $p.value
# [1] 1.113095e-20

c_index_mat <- sapply(levels(clin_comp$clust), function(x){
    clin = concordance(Surv(OS_YEARS, OS_STATUS) ~ clin_score, data = clin_comp, 
        subset = clust == x )$concordance
    mol = concordance(Surv(OS_YEARS, OS_STATUS) ~ mol_score, data = clin_comp, 
        subset = clust == x )$concordance
    c(IPSSM_clinical = clin, IPSSM = mol)
})

png("figures/tsneclust_cindex.png")
c_index_mat %>% 
    t() %>%
    as_tibble() %>%
    mutate(cluster = factor(1:17)) %>%
    pivot_longer(cols = c(IPSSM_clinical, IPSSM)) %>%
    ggplot(aes(x = cluster, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw()
dev.off()




## Interacciones
#' ---------------------------
#' 
#' ## Sex
#' ---------------------------
cox_sex_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, clinical_clust, subset = clust != 14)
coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*hma + AGE + SEX, clinical_clust)

sex_int_clust <- subset(sex_prop_df, propF > 0.25 & propF < 0.75)$cluster

cox_sex_int <- clinical_clust %>%
    subset(clust %in% sex_int_clust) %>% ## Remove clusters with high sex-imbalance (>1:3)
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SEX + AGE, data = .)



sex_int_coef <- cox_sex_int$coef
sex_int_coef <- sex_int_coef[grep(":SEXM", names(sex_int_coef))]
sex_int_coef <- c("clust1:SEXM" = 0, sex_int_coef)
sex_int_coef <- gsub("clust", "", names(sort(sex_int_coef)))
sex_int_coef <- gsub(":SEXM", "", sex_int_coef)

png("figures/tsne_surv_clustSex.png", height = 800, width = 1200)
a <- lapply(sex_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clinical_clust, subset = clust == x) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    subset(clust %in% sex_int_clust) %>% ## Remove clusters with high sex-imbalance
    mutate(superClust = ifelse(clust %in% c(10, 3, 1), "equal", "male")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*superClust + AGE, data = .)

clinical_clust %>%
    subset(clust %in% c("3", "5")) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*clust + AGE, data = .)


lapply(sex_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX + AGE, 
        clinical_clust, subset = clust == x) 
})

#' ## IPSSM
#' ---------------------------

cox_ipssm_int <- coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*IPSSM_SCORE + SEX + AGE, data = clinical_clust)
cox_ipssra_int <-coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*IPSSRA_SCORE + SEX + AGE, data = clinical_clust)


cox_ipssm_int_coef <- cox_ipssm_int$coef
cox_ipssra_int_coef <- cox_ipssra_int$coef

cox_ipss_int_df <- tibble(cluster = factor(paste0("clust", c(7, 1:6, 8:17)), levels = paste0("clust", 1:17)),
    IPSSM = c(0, cox_ipssm_int_coef[grep(":IPSSM_SCORE", names(cox_ipssm_int_coef))]),
    IPSSRA = c(0, cox_ipssra_int_coef[grep(":IPSSRA_SCORE", names(cox_ipssra_int_coef))]))

png("figures/tsne_clust_ipss_int.png", width = 800)
cox_ipss_int_df %>%
    pivot_longer(!cluster, names_to = "Score", values_to = "Coef") %>%
    mutate(Score = factor(Score, levels = c("IPSSRA", "IPSSM"))) %>%
    ggplot(aes(x = cluster, y = Coef, fill = Score)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw()
dev.off()

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
    filter(clust %in% c(8, 5, 3, 4, 10, 11, 13, 2)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

hma_int_coef_h <- cox_hma_int_h$coef
hma_int_coef_h <- hma_int_coef_h[grep("hma:", names(hma_int_coef_h))]
hma_int_coef_h <- c("hma:clust2" = 0, hma_int_coef_h)
hma_int_coef_h <- gsub("hma:clust", "", names(sort(hma_int_coef_h)))
png("figures/tsne_surv_clusthma_high.png", width = 800)
a <- lapply(hma_int_coef_h, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    ggsurvplot(data = clinical_clust, title = paste0("Cluster", x))
    p$plot
})
plot_grid(plotlist = a)
dev.off()

clinical_clust %>%
    filter(clust %in% c(8, 5, 3, 4, 10, 11, 13, 2)) %>%
    mutate(superClust = ifelse(clust %in% c(2, 3, 10, 8), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*superClust + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

names(hma_int_coef_h) <- hma_int_coef_h
lapply(hma_int_coef_h, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE, 
        clinical_clust, subset = clust == x & IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) 
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
    filter(clust %in% c(6, 5, 3, 8, 11, 13)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + SEX + AGE, data = ., subset = IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)
hma_int_coef_vh <- cox_hma_int_vh$coef
hma_int_coef_vh <- hma_int_coef_vh[grep("hma:", names(hma_int_coef_vh))]
hma_int_coef_vh <- c("hma:clust3" = 0, hma_int_coef_vh)
hma_int_coef_vh <- gsub("hma:clust", "", names(sort(hma_int_coef_vh)))
png("figures/tsne_surv_clusthma_veryhigh.png", width = 800)
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


png("figures/tsne_surv_clusthma_veryhigh_clustGroup.png")
cot <- clinical_clust %>%
     subset(IPSSM == "Very-High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0) %>%
    mutate(sC = ifelse(clust %in% c(6, 3), "low", "high")) 
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
    subset(clust %in% c(5, 11, 3, 6, 8) & IPSSM %in% c("High", "Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6) & IPSSM %in% c("High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


clinical_clust %>%
   subset(clust %in% c(5, 11, 3, 6, 8) & IPSSM %in% c("Very-High")) %>%
    mutate(sc = ifelse(clust %in% c(3, 6), "low", "high")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*sc + SEX + AGE, data = ., subset =  lenalidomid == 0 & chemotherapy == 0 & transplant == 0)

#' Parece que el cluster 3 es menos efectivo a hma
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

                                              "MDS/MPN-U","other"))
dd %>%
   subset(!who %in% c("other", "MDS/MPN-RS-T", "MDS/MPN-U")) %>%
   mutate(who = droplevels(who)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*who + SEX + AGE, data = ., subset = IPSSM == "High" & lenalidomid == 0 & chemotherapy == 0 & transplant == 0)


#' ---------------------------
#' # Mutations
#' ---------------------------


#' # SF3B1
#' ---------------------------


sf3b1_props <- prop.table(table(clin_mut$SF3B1, clin_mut$clust), margin = 2)[2, ]
sf3b1_clust <- names(sf3b1_props)[sf3b1_props > 0.1] ## Select clusters with > 10% of SF3B1 

cox_sf3b1_int <- clin_mut %>%
    subset(clust %in% sf3b1_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SF3B1 + SEX + AGE + IPSSR_SCORE, data = .)


sf3b1_int_coef <- cox_sf3b1_int$coef
sf3b1_int_coef <- sf3b1_int_coef[grep(":SF3B1", names(sf3b1_int_coef))]
sf3b1_int_coef <- c("clust2:SF3B1" = 0, sf3b1_int_coef)
sf3b1_int_coef <- gsub("clust", "", names(sort(sf3b1_int_coef)))
sf3b1_int_coef <- gsub(":SF3B1", "", sf3b1_int_coef)

png("figures/tsne_surv_clustSF3B1.png", height = 800, width = 1200)
a <- lapply(sf3b1_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SF3B1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

lapply(sf3b1_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SF3B1  + AGE, 
        clin_mut, subset = clust == x) 
})


#' # ASXL1
#' ---------------------------

asxl1_props <- prop.table(table(clin_mut$ASXL1, clin_mut$clust), margin = 2)[2, ]
asxl1_clust <- names(asxl1_props)[asxl1_props > 0.1] ## Select clusters with > 10% of ASXL1 

cox_asxl1_int <- clin_mut %>%
    subset(clust %in% asxl1_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*ASXL1 + SEX + AGE, data = .)


asxl1_int_coef <- cox_asxl1_int$coef
asxl1_int_coef <- asxl1_int_coef[grep(":ASXL1", names(asxl1_int_coef))]
asxl1_int_coef <- c("clust1:ASXL1" = 0, asxl1_int_coef)
asxl1_int_coef <- gsub("clust", "", names(sort(asxl1_int_coef)))
asxl1_int_coef <- gsub(":ASXL1", "", asxl1_int_coef)

png("figures/tsne_surv_clustASXL1.png", height = 800, width = 1200)
a <- lapply(asxl1_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    if (x != "8") {
      p$plot + theme_bw(base_size = 15)  + theme(legend.position = "none") 
    } else {
        p$plot + theme_bw(base_size = 15) 
    }
})
plot_grid(plotlist = a)
dev.off()

lapply(asxl1_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1  + AGE, 
        clin_mut, subset = clust == x) 
})


#' # DNMT3A
#' ---------------------------
dnmt3a_props <- prop.table(table(clin_mut$DNMT3A, clin_mut$clust), margin = 2)[2, ]
dnmt3a_clust <- names(dnmt3a_props)[dnmt3a_props > 0.10] ## Select clusters with > 15% of DNMT3A 

cox_dnmt3a_int <- clin_mut %>%
    subset(clust %in% dnmt3a_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*DNMT3A + SEX + AGE, data = .)

dnmt3a_int_coef <- cox_dnmt3a_int$coef
dnmt3a_int_coef <- dnmt3a_int_coef[grep(":DNMT3A", names(dnmt3a_int_coef))]
dnmt3a_int_coef <- c("clust1:DNMT3A" = 0, dnmt3a_int_coef)
dnmt3a_int_coef <- gsub("clust", "", names(sort(dnmt3a_int_coef)))
dnmt3a_int_coef <- gsub(":DNMT3A", "", dnmt3a_int_coef)

png("figures/tsne_surv_clustDNMT3A.png", height = 800, width = 1200)
a <- lapply(dnmt3a_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

lapply(dnmt3a_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SF3B1  + AGE, 
        clin_mut, subset = clust == x) 
})

#' # SRSF2
#' ---------------------------
srsf2_props <- prop.table(table(clin_mut$SRSF2, clin_mut$clust), margin = 2)[2, ]
srsf2_clust <- names(srsf2_props)[srsf2_props > 0.1] ## Select clusters with > 10% of SRSF2 

cox_srsf2_int <- clin_mut %>%
    subset(clust %in% srsf2_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*SRSF2 + SEX + AGE, data = .)


srsf2_int_coef <- cox_srsf2_int$coef
srsf2_int_coef <- srsf2_int_coef[grep(":SRSF2", names(srsf2_int_coef))]
srsf2_int_coef <- c("clust1:SRSF2" = 0, srsf2_int_coef)
srsf2_int_coef <- gsub("clust", "", names(sort(srsf2_int_coef)))
srsf2_int_coef <- gsub(":SRSF2", "", srsf2_int_coef)

png("figures/tsne_surv_clustSRSF2.png", height = 800, width = 1200)
a <- lapply(srsf2_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SRSF2, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

lapply(srsf2_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SF3B1  + AGE, 
        clin_mut, subset = clust == x) 
})


#' # RUNX1
#' ---------------------------
runx1_props <- prop.table(table(clin_mut$RUNX1, clin_mut$clust), margin = 2)[2, ]
runx1_clust <- names(runx1_props)[runx1_props > 0.1] ## Select clusters with > 10% of RUNX1 

cox_runx1_int <- clin_mut %>%
    subset(clust %in% runx1_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*RUNX1 + SEX + AGE, data = .)


runx1_int_coef <- cox_runx1_int$coef
runx1_int_coef <- runx1_int_coef[grep(":RUNX1", names(runx1_int_coef))]
runx1_int_coef <- c("clust1:RUNX1" = 0, runx1_int_coef)
runx1_int_coef <- gsub("clust", "", names(sort(runx1_int_coef)))
runx1_int_coef <- gsub(":RUNX1", "", runx1_int_coef)

png("figures/tsne_surv_clustRUNX1.png", height = 800, width = 1200)
a <- lapply(runx1_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

lapply(runx1_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1  + AGE, 
        clin_mut, subset = clust == x) 
})


#' # BCOR
#' ---------------------------
bcor_props <- prop.table(table(clin_mut$BCOR, clin_mut$clust), margin = 2)[2, ]
bcor_clust <- names(bcor_props)[bcor_props > 0.05] ## Select clusters with > 10% of RUNX1 

cox_bcor_int <- clin_mut %>%
    subset(clust %in% bcor_clust) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*BCOR + SEX + AGE, data = .)


bcor_int_coef <- cox_bcor_int$coef
bcor_int_coef <- bcor_int_coef[grep(":BCOR", names(bcor_int_coef))]
bcor_int_coef <- c("clust3:BCOR" = 0, bcor_int_coef)
bcor_int_coef <- gsub("clust", "", names(sort(bcor_int_coef)))
bcor_int_coef <- gsub(":BCOR", "", bcor_int_coef)

png("figures/tsne_surv_clustBCOR.png", height = 800, width = 1200)
a <- lapply(bcor_int_coef, function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ BCOR, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

lapply(bcor_int_coef, function(x){
    print(x)
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1  + AGE, 
        clin_mut, subset = clust == x) 
})




## Report 2024
#' ---------------------------

kar_plot <- pheatmap(mat = t(kar_tabs_prop2[, c("complex", "monoy", "del5q", "mono7", "triso8")]), 
    display_numbers = TRUE, font_size = 20 )

png("figures/eval_2024/cluster_karyotipe_mini.png", height = 900, width = 2100, res = 300)
kar_plot
dev.off()

surv_plot <- ggsurvplot(os_clust, data = clinical_clust, 
    palette = c("grey", "grey", "coral", "grey", "maroon", "red", "cyan", 
        "pink", rep("grey", 5), "blue", rep("grey", 2),  "aquamarine"))

png("figures/eval_2024/cluster_survival.png", height = 1800, width = 2100, res = 300)
surv_plot
dev.off()



who_types <- pheatmap(t(prop.table(clust_who_tab, margin = 1)*100), display_numbers = TRUE)

png("figures/eval_2024/cluster_who_mini.png", height = 1300, width = 2400, res = 300)
who_types
dev.off()


mut_plot <- pheatmap(t(mut_tabs_prop2[, c("mutation0", "mutation1", "mutation2_3", 
    "mutation_multi", "TP53multi", "SF3B1", "RUNX1", "NRAS", "ETV6",
    "IDH2", "CBL", "EZH2", "U2AF1", "SRSF2", "DNMT3A", "ASXL1", "KRAS", "TET2")]), display_numbers = TRUE)

png("figures/eval_2024/cluster_mutation_mini.png", height = 2000, width = 2400, res = 300)
mut_plot
dev.off()


png("figures/eval_2024/features1_panel.png", width = 1000, height = 1000)
plot_grid(
    plot_grid(age_plot + theme_bw(base_size = 20),kar_plot$gtable , nrow = 1, labels = "AUTO"), 
    plot_grid(surv_plot$plot, who_types$gtable, nrow = 1, rel_widths = c(1, 1.5), labels = c("C", "D")),
    mut_plot$gtable, 
nrow = 3, rel_heights = c(1, 1.4, 1.7), labels = c("", "", "E"))
dev.off()


## Interactions with sex
sex_hr_clust <- lapply(1:17, function(x){
        freq <- mean(clin_mut$SEX[clin_mut$clust == x] == "M")
        if (freq > 0.05 & freq < 0.95){
            mod <- coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX + AGE + IPSSM_SCORE, 
             clin_mut, subset = clust == x)
            coef <- summary(mod)$coefficients
            hr <- exp(coef[1, 1])
            pval <- coef[1, 5]
        } else{
            hr <- 1
            pval <- 1
        }
        data.frame(HR = hr, Pvalue = pval, Cluster = x, Freq = freq )
    }) %>% Reduce(f = rbind) %>%
        as_tibble() %>%
    mutate(Cluster = factor(Cluster, levels = 1:17))


sex_hr_plot <-  mutate(sex_hr_clust, 
    Sig = ifelse(Pvalue < 0.05, "Signif", "Not-signif"),
    Direction = ifelse(HR > 1, "Risk", "Protective"), 
    Color = ifelse(Freq < 0.1 | Freq > 0.9, "Undetermined", paste(Direction, Sig)),
    Color = factor(Color, levels = c("Protective Signif", "Protective Not-signif",
        "Risk Signif", "Risk Not-signif", "Undetermined")) ) %>%
    filter(Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey")) +
        ggtitle("HR for being male") +
        theme(plot.title = element_text(hjust = 0.5))

png("figures/eval_2024/cluster_interaction_sex.png", height = 1000, width = 1800, res = 300)
sex_hr_plot
dev.off()


clin_mut %>%
    filter(clust %in% c(2, 9, 10, 16)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*clust + AGE + IPSSM_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(9, 10, 16)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*clust + AGE + IPSSM_SCORE, data = .) 



clin_sex_clust <- clin_mut %>%
    filter(!clust %in% c(4, 11, 14)) %>%
    mutate(clust = droplevels(clust),
            sex_clust = ifelse(clust %in% c(7, 10, 16), "Female", 
                ifelse(clust %in% c(2, 9), "Male", "Baseline")),
            sex_clust = factor(sex_clust, levels = c("Baseline", "Male","Female")))

coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX*sex_clust + AGE + IPSSM_SCORE, data = clin_sex_clust)

lapply(c("Female", "Male", "Baseline"), function(x){
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX + AGE + IPSSM_SCORE, clin_sex_clust, subset = sex_clust == x)
})

coxph(Surv(OS_YEARS,OS_STATUS) ~ SEX + AGE + IPSSM_SCORE, clin_sex_clust)

png("figures/eval_2024/cluster_interaction_sexgroup_surv.png", width = 4000, height = 1000, res = 300)
a <- lapply(c("Female", "Baseline", "Male"), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clin_sex_clust, subset = sex_clust == x) %>%
    ggsurvplot(data = clin_mut, title = x)
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a, nrow = 1)
dev.off()


png("figures/eval_2024/cluster_interaction_sex_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(2, 9, 10, 16), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SEX, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()



## Interactions with scores
clinical_clust$IPSSMOL_score <- clinical_clust$IPSSM_SCORE - clinical_clust$IPSSMscore
scores_hr_clust <-  lapply(1:17, function(x){
        mod <- coxph(Surv(OS_YEARS,OS_STATUS) ~ IPSSMscore + IPSSMOL_score + AGE, 
            clinical_clust, subset = clust == x)
        coef <- summary(mod)$coefficients
        hrs <- exp(coef[1:2, 1])
        pvals <- coef[1:2, 5]
        data.frame(HR = hrs, Pvalue = pvals, Cluster = x, Score = c("Clinical", "Molecular") )
    }) %>% Reduce(f = rbind) %>%
        as_tibble() %>%
    mutate(Cluster = factor(Cluster, levels = 1:17))

scores_hr_plot <-  scores_hr_clust %>%
    ggplot(aes(x = Cluster, y = HR, fill = Score)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() +
        scale_y_log10() +
        ggtitle("IPSSM decomposition") +
        theme(plot.title = element_text(hjust = 0.5))


## Interactions with treatment
treat_df <- filter(clinical_clust, lenalidomid == 0 & chemotherapy == 0 & transplant == 0 & !is.na(IPSSM))

hma_hr_clust <- lapply(c("High", "Very-High"), function(risk){
    lapply(1:17, function(x){
        N <- sum(treat_df$clust == x & treat_df$IPSSM == risk)
        freq <- mean(treat_df$hma[treat_df$clust == x & treat_df$IPSSM == risk])
        if (freq == 0){
            hr <- 1
            pval <- 1
        } else{
            mod <- coxph(Surv(OS_YEARS,OS_STATUS) ~ hma + AGE + IPSSM_SCORE, 
                treat_df, subset = clust == x & IPSSM == risk)
            coef <- summary(mod)$coefficients
            hr <- exp(coef[1, 1])
            pval <- coef[1, 5]
        }
        data.frame(HR = hr, Pvalue = pval, Cluster = x, Risk = risk, Freq = freq,
        N = N )
    }) %>% Reduce(f = rbind) %>%
        as_tibble()
}) %>% Reduce(f = rbind)

hma_hr_col <- mutate(hma_hr_clust, 
    Sig = ifelse(Pvalue < 0.05, "Signif", "Not-signif"),
    Direction = ifelse(HR > 1, "Risk", "Protective"), 
    Color = ifelse(Freq < 0.1 | N < 30, "Undetermined", paste(Direction, Sig)),
    Color = factor(Color, levels = c("Protective Signif", "Protective Not-signif",
        "Risk Signif", "Risk Not-signif", "Undetermined")) )


sel_clusters <- group_by(hma_hr_col, Cluster) %>%
    summarize(N = sum(Color != "Undetermined")) %>%
    filter(N > 0) %>%
    pull(Cluster)

hma_plot <- hma_hr_col %>%
    filter(Cluster %in% sel_clusters) %>%
    mutate(Cluster = factor(Cluster, levels = sel_clusters),
        HR = ifelse(Color == "Undetermined", 1, HR)) %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        facet_wrap(~ Risk) +
        scale_fill_manual(name = "", values = c("green", "darkgreen", "red", "darkred", "grey")) +
        ggtitle("HMA treatment") +
        theme(plot.title = element_text(hjust = 0.5))



png("figures/eval_2024/hma_cluster_interaction.png", height = 1000, width = 1800, res = 300)
hma_plot
dev.off()

png("figures/eval_2024/hma_high_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(3, 8, 5, 4), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clin_mut, subset = clust == x & IPSSM == "High") %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()


png("figures/eval_2024/hma_veryhigh_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(3, 6, 5, 8), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ hma, clin_mut, subset = clust == x & IPSSM == "Very-High") %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(3, 4, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "4")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("High", "Very-High")) 

clin_mut %>%
    filter(clust %in% c(3, 4, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "5")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("High")) 

clin_mut %>%
    filter(clust %in% c(3, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "3")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("Very-High")) 


clin_mut %>%
    filter(clust %in% c(3, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "3")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("Very-High")) 

clin_mut %>%
    filter(clust %in% c(3, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "6")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("Very-High")) 


clin_mut %>%
    filter(clust %in% c(3, 6, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "5")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ hma*clust + AGE + SEX + IPSSM_SCORE, data = ., subset = IPSSM %in% c("Very-High")) 



## Interactions with mutations
sel_genes <- c("RUNX1", "SRSF2", "DNMT3A", "ASXL1", "SF3B1", "TET2")
mut_hr_clust <- lapply(sel_genes, function(gene){
    lapply(1:17, function(x){
        mod <- coxph(formula(paste("Surv(OS_YEARS,OS_STATUS) ~", gene, " + AGE + IPSSR_SCORE")), 
            clin_mut, subset = clust == x)
        coef <- summary(mod)$coefficients
        hr <- exp(coef[1, 1])
        pval <- coef[1, 5]
        freq <- mean(clin_mut[[gene]][clin_mut$clust == x])
        data.frame(HR = hr, Pvalue = pval, Cluster = x, Gene = gene, Freq = freq )
    }) %>% Reduce(f = rbind) %>%
        as_tibble()
}) %>% Reduce(f = rbind) %>%
    mutate(Cluster = factor(Cluster, levels = 1:17))

## Modify one value due to low frequency
mut_hr_clust <- mut_hr_clust %>%
    mutate(HR = ifelse(HR < 0.41, 0.41, HR))

## add colors
mut_hr_col <- mutate(mut_hr_clust, 
    Sig = ifelse(Pvalue < 0.05, "Signif", "Not-signif"),
    Direction = ifelse(HR > 1, "Risk", "Protective"), 
    Color = ifelse(Freq < 0.1, "Undetermined", paste(Direction, Sig)),
    Color = factor(Color, levels = c("Protective Signif", "Protective Not-signif",
        "Risk Signif", "Risk Not-signif", "Undetermined")) )


mut_hr_plot <- mut_hr_col %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        facet_wrap(~ Gene) +
        scale_fill_manual(name = "", values = c("green", "darkgreen", "red", "darkred", "grey"))


png("figures/eval_2024/mutations_risk.png", width = 3000, height = 1000, res = 300)
mut_hr_plot
dev.off()


png("figures/eval_2024/interactions.png", width = 1000, height = 500)
plot_grid(
    plot_grid(scores_hr_plot, sex_hr_plot, hma_plot, ncol = 3, labels = c("A", "B", "D")),
    mut_hr_plot,
nrow = 2, labels = c("", "C"))

dev.off()


#' # SF3B1
#' ---------------------------
sf3b1_plot <- mut_hr_col %>%
    filter(Gene == "SF3B1" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey"))
png("figures/eval_2024/sf3b1_risk.png", width = 2000, height = 800, res = 300)
sf3b1_plot
dev.off()

png("figures/eval_2024/sf3b1_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(15, 9, 12, 17), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SF3B1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(15, 9, 12, 17)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SF3B1*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(15, 9, 12, 17)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "15")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SF3B1*clust + AGE + IPSSR_SCORE, data = .) 



#' # ASXL1
#' ---------------------------
asxl1_plot <- mut_hr_col %>%
    filter(Gene == "ASXL1" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey"))
png("figures/eval_2024/asxl1_risk.png", width = 2000, height = 800, res = 300)
asxl1_plot
dev.off()

png("figures/eval_2024/asxl1_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(17, 7, 9, 16), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ ASXL1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(17, 7, 9, 16)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(17, 7, 9, 16)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "17")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ ASXL1*clust + AGE + IPSSR_SCORE, data = .) 



#' # DNMT3A
#' ---------------------------
dnmt3a_plot <- mut_hr_col %>%
    filter(Gene == "DNMT3A" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("green", "darkgreen", "darkred", "grey"))
png("figures/eval_2024/dnmt3a_risk.png", width = 2000, height = 800, res = 300)
dnmt3a_plot
dev.off()

png("figures/eval_2024/dnmt3a_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(1, 3, 11, 13), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ DNMT3A, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(1, 3, 11, 13)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(3, 11, 13)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ DNMT3A*clust + AGE + IPSSR_SCORE, data = .) 



#' # RUNX1
#' ---------------------------
runx1_plot <- mut_hr_col %>%
    filter(Gene == "RUNX1" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("red", "darkred", "grey"))
png("figures/eval_2024/runx1_risk.png", width = 2000, height = 800, res = 300)
runx1_plot
dev.off()

png("figures/eval_2024/runx1_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(11, 10, 3, 8), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ RUNX1, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(11, 10, 3, 8)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(11, 10, 8)) %>%
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ RUNX1*clust + AGE + IPSSR_SCORE, data = .) 



#' # SRSF2
#' ---------------------------
srsf2_plot <- mut_hr_col %>%
    filter(Gene == "SRSF2" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey"))
png("figures/eval_2024/srsf2_risk.png", width = 2000, height = 800, res = 300)
srsf2_plot
dev.off()

png("figures/eval_2024/srsf2_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(4, 7, 1, 13), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ SRSF2, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(4, 7, 1, 13)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "4")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(4, 7, 1, 13)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "7")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ SRSF2*clust + AGE + IPSSR_SCORE, data = .) 




#' # TET2
#' ---------------------------
tet2_plot <- mut_hr_col %>%
    filter(Gene == "TET2" & Color != "Undetermined") %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey"))
png("figures/eval_2024/tet2_risk.png", width = 2000, height = 800, res = 300)
tet2_plot
dev.off()

png("figures/eval_2024/tet2_cluster_surv.png", width = 3000, height = 2000, res = 300)
a <- lapply(c(17, 12, 1, 9), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ TET2, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a)
dev.off()

clin_mut %>%
    filter(clust %in% c(17, 12, 1, 9)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "17")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ TET2*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(17, 12, 1, 9)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "12")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ TET2*clust + AGE + IPSSR_SCORE, data = .) 


#' # BCOR
#' ---------------------------
#' 
#' 

bcor_hr_clust <- lapply(1:17, function(x){
        mod <- coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR + AGE + IPSSR_SCORE, 
            clin_mut, subset = clust == x)
        coef <- summary(mod)$coefficients
        hr <- exp(coef[1, 1])
        pval <- coef[1, 5]
        freq <- mean(clin_mut[["BCOR"]][clin_mut$clust == x])
        data.frame(HR = hr, Pvalue = pval, Cluster = x, Gene = "BCOR", Freq = freq )
    }) %>% Reduce(f = rbind) %>%
        as_tibble() %>%
    mutate(Cluster = factor(Cluster, levels = 1:17),
    Sig = ifelse(Pvalue < 0.05, "Signif", "Not-signif"),
    Direction = ifelse(HR > 1, "Risk", "Protective"), 
    Color = ifelse(Freq < 0.05, "Undetermined", paste(Direction, Sig)),
    Color = factor(Color, levels = c("Protective Signif", "Protective Not-signif",
        "Risk Signif", "Risk Not-signif", "Undetermined")) )



bcor_plot <- bcor_hr_clust %>%
    filter(Freq > 0.05) %>%
    ggplot(aes(x = Cluster, y = HR, fill = Color)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        scale_fill_manual(name = "", values = c("darkgreen", "red", "darkred", "grey"))
png("figures/eval_2024/bcor_risk.png", width = 2000, height = 800, res = 300)
bcor_plot
dev.off()

png("figures/eval_2024/bcor_cluster_surv.png", width = 3000, height = 2500, res = 300)
a <- lapply(c(4, 11, 3, 5, 8), function(x){
p <-  survfit(formula = Surv(OS_YEARS,OS_STATUS) ~ BCOR, clin_mut, subset = clust == x) %>%
    ggsurvplot(data = clin_mut, title = paste0("Cluster", x))
    p$plot + theme_bw(base_size = 20) 
})
plot_grid(plotlist = a, ncol = 2)
dev.off()

clin_mut %>%
    filter(clust %in% c(4, 11, 3, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "4")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(4, 11, 3, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "11")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR*clust + AGE + IPSSR_SCORE, data = .) 

clin_mut %>%
    filter(clust %in% c(4, 11, 3, 5, 8)) %>%
    mutate(clust = droplevels(clust),
            clust = relevel(clust, "3")) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR*clust + AGE + IPSSR_SCORE, data = .) 




clin_mut %>%
    subset(clust %in% c(3, 5, 8)) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*BCOR + SEX + AGE + IPSSR_SCORE, data = .)

coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR + SEX + AGE + IPSSR_SCORE, data = clin_mut, subset = clust == 3)
coxph(Surv(OS_YEARS,OS_STATUS) ~ BCOR + SEX + AGE + IPSSR_SCORE, data = clin_mut, subset = clust == 5)






clin_mut %>%
    subset(clust %in% c(3, 5)) %>% 
    mutate(clust = droplevels(clust)) %>%
    coxph(Surv(OS_YEARS,OS_STATUS) ~ clust*U2AF1 + SEX + AGE + IPSSR_SCORE, data = .)

coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + SEX + AGE + IPSSR_SCORE, data = clin_mut, subset = clust == 3)
coxph(Surv(OS_YEARS,OS_STATUS) ~ U2AF1 + SEX + AGE + IPSSR_SCORE, data = clin_mut, subset = clust == 5)

