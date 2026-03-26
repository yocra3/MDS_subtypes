#' ---------------------------
#'
#' Purpose of script:
#'
#'  Run transcriptomic analysis using limma between sub-groups
#' 
#' ---------------------------
#'
#' Notes:
#' Make figures for subgroups exploration
#' 
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.8 R
#'
#' ---------------------------


## Load libraries
library(limma)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)
library(NetActivity)
library(pheatmap)

## Load data
load("results/hershberger/rnaseq_hersh.Rdata")

## Explore raw data
vst_hersh <- vst(DESeqDataSet(se_hersh, design = ~ WHO), blind = TRUE)

vst_hersh_nosex <- vst_hersh[!seqnames(vst_hersh) %in% c("X", "Y"), ]

pdf("figures/rnaseq_hersh/pca_raw.pdf")
plotPCA(vst_hersh_nosex, intgroup = "WHO") + theme_bw()
plotPCA(vst_hersh_nosex, intgroup = "SEX") + theme_bw()
dev.off()

## Select only MDS samples
se_hersh_mds <- se_hersh[, grepl("MDS", se_hersh$WHO) ]
vst_hersh_mds <- vst(DESeqDataSet(se_hersh_mds, design = ~ WHO), blind = FALSE)
vst_hersh_mds_nosex <- vst_hersh_mds[!seqnames(vst_hersh_mds) %in% c("X", "Y"), ]


pdf("figures/rnaseq_hersh/pca_mds.pdf")
plotPCA(vst_hersh_mds_nosex, intgroup = "WHO") + theme_bw()
plotPCA(vst_hersh_mds_nosex, intgroup = "SEX") + theme_bw()
plotPCA(vst_hersh_mds_nosex, intgroup = "sub_group") + theme_bw()
dev.off()


png("figures/rnaseq_hersh/pca_mds_subgroup.png")
plotPCA(vst_hersh_mds_nosex, intgroup = "sub_group") + theme_bw() + 
scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", 
    "#D55E00", "#999999", "grey40",  "black"))
dev.off()

## Pot distance heatmap
sampleDists <- dist(t(assay(vst_hersh_mds_nosex)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

col_annot <- data.frame(sub_group = vst_hersh_mds_nosex$sub_group)
rownames(col_annot) <- colnames(vst_hersh_mds_nosex)

png("figures/rnaseq_hersh/heatmap_distance.png", width = 800, height = 800)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         annotation_col = col_annot)
dev.off()

clust <- cutree(hclust(sampleDists), k = 20)

## Filter by gene expression
se_hersh_mds$sub_group <- relevel(factor(se_hersh_mds$sub_group), ref = "Low blasts")
design <- model.matrix(~ 1 + sub_group, colData(se_hersh_mds))

keep <- filterByExpr(se_hersh_mds, design)

se_hersh_filt <- se_hersh_mds[keep, ]


## Transform with voom
se_hersh_filt <- calcNormFactors(se_hersh_filt)
se_hersh_voom <- voom(se_hersh_filt, design)

## Fit linear model
lm_design <- model.matrix(~ 1 + sub_group + AGE + SEX, colData(se_hersh_mds))

main_fit <- lmFit(se_hersh_voom, lm_design)
main_fit <- eBayes(main_fit)

tables <- lapply(2:8, function(i){
    topTable(main_fit, coef = i, n = Inf)
})
names(tables) <- levels(se_hersh_mds$sub_group)[-c(1, 9:10)]


lapply(tables, function(x) filter(x, adj.P.Val < 0.05))
sapply(tables, function(x) sum(x$adj.P.Val < 0.05))

## NetActivity
vst_hersh_mds_ens <- vst_hersh_mds
rownames(vst_hersh_mds_ens) <- rowRanges(vst_hersh_mds_ens)$gene_id
net_hersh_mds <- prepareSummarizedExperiment(vst_hersh_mds_ens, "gtex_gokegg")
net_hersh_mds <- net_hersh_mds[!duplicated(rownames(net_hersh_mds)), ]
net_hersh_mds_scores <- computeGeneSetScores(net_hersh_mds, "gtex_gokegg")


net_fit <- lmFit(assay(net_hersh_mds_scores), lm_design)
net_fit <- eBayes(net_fit)

net_tables <- lapply(2:8, function(i){
    topTable(net_fit, coef = i, n = Inf)
})
names(net_tables) <- levels(se_hersh_mds$sub_group)[-c(1, 9:10)]
net_tables_filt <- lapply(net_tables, function(x) filter(x, adj.P.Val < 0.05))

net_tables_filt <- lapply(net_tables, function(x) filter(x, P.Value < 1e-6))

sapply(net_tables, function(x) sum(x$adj.P.Val < 0.05))
