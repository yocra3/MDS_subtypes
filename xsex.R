library(GEOquery)
library(tidyverse)
library(minfi)
library(readxl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

geo <- getGEO("GSE152710")

mat <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152710/suppl/GSE152710%5FMatrix%5Fprocessed.tsv.gz")

mat_filt <- select(mat, -ends_with("pval"))
mat_cpg <- data.matrix(mat_filt[, -1])
rownames(mat_cpg) <- mat_filt$ID_REF

map <- pData(geo[[1]])[, c("geo_accession", "description")]
rownames(map) <- map$description

colnames(mat_cpg) <- map[colnames(mat_cpg), "geo_accession"]

gset <- makeGenomicRatioSetFromMatrix(mat_cpg[, map$geo_accession], pData = pData(geo[[1]]))
male_inds <- mat[mat$ID_REF == "cg09595415",] %>% select(ends_with("pval")) %>% unlist() < 2e-16
male_ids <- colnames(mat_cpg)[male_inds]
gset$sex <- ifelse(gset$geo_accession %in% male_ids, "male", "female")

pseudo_reg <- GRanges(c("chrX:60001-2699520", "chrX:154931044-155260560"))
chrX <- subset(rowRanges(gset), seqnames == "chrX")
chrX_filt <- subsetByOverlaps(chrX, pseudo_reg, invert = TRUE)
mat[mat$ID_REF == "cg08265308",] %>% select(ends_with("pval"))

df <- data.frame(beta = getBeta(gset["cg17939569",]) %>% as.numeric(), sex = gset$sex)
boxplot(beta ~ sex, df)
dev.off()

gset_female <- gset[, gset$sex == "female"]

xci_genes <- read_xlsx("data/XCI_genes.PMID29022598.xlsx", skip = 1) 
#genes that always escape
always <- xci_genes$`Gene name`[xci_genes$`Combined XCI status`=="escape"]
#genes that are inactive
inactive <- xci_genes$`Gene name`[xci_genes$`Combined XCI status`=="inactive"]

gset_female_chrX <- subsetByOverlaps(gset_female, GRanges("chrX:1-999999999"))
data(Other)
rowData(gset_female_chrX) <- Other[rownames(gset_female_chrX), ]
gset_inact <- gset_female_chrX[rowData(gset_female_chrX)$UCSC_RefGene_Name %in% inactive,]

RXA <- colMeans(getBeta(gset_inact) < 0.2)

df <- data.frame(RXA = RXA, response = factor(gset_inact$`response:ch1`), 
treatment = gset_female$`treatment:ch1`)
glm(response ~ RXA + treatment, df, family = "binomial", subset = response != "NA") %>% summary()
boxplot(RXA ~ response, df)
dev.off()


geo2 <- getGEO("GSE124413")

mat2 <- read_csv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124413/suppl/GSE124413%5Fmatrix%5Fprocessed.csv.gz")