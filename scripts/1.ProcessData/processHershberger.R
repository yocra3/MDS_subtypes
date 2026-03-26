#' ---------------------------
#'
#' Purpose of script:
#'
#'  Process data from Hershberger to compare with our results
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.8 R
#'
#' ---------------------------
#' 

# Load libraries and data
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(SummarizedExperiment)
library(readxl)
library(ipssm)

## Load data
samples <- read_xlsx("data/hershberger/20250904_Sample_Information_MDS_865.xlsx")
os <- read_xlsx("data/hershberger/20251022_Survival_Information_MDS_865.xlsx")
mutations <- read_xlsx("data/hershberger/20250905_Mutations_MDS_865.xlsx") 
tet2bi <- read_xlsx("data/hershberger/20260305_TET2_MDS_share.xlsx") 
tet2bi_v1 <- read_xlsx("data/hershberger/20260223_TET2_MDS_share.xlsx") 

tet2bi <- rbind(tet2bi, tet2bi_v1)
check_cols <- function(event, event2, var1, var2, var3, var4, var5){


    boolean <- grepl(event, var1, fixed = TRUE) 
    boolean2 <- grepl(event2, var2, fixed = TRUE) | 
        grepl(event2, var3, fixed = TRUE) | grepl(event2, var4, fixed = TRUE) | grepl(event2, var5, fixed = TRUE)
    ifelse(boolean | boolean2, 1, 0)


}

hershberger <- samples %>%
    mutate(ID = `Exam Array`,
    WHO = `WHO 2022 Diagnosis`, 
    AGE = Age, 
    SEX = ifelse(Sex == "male", "M", "F"),
    BM_BLAST = as.numeric(`% Blasts im BM`), 
    WBC = as.numeric(`WBC in µl`)/1000, HB = as.numeric(`HB in g/dl`), 
    PLT = as.numeric(`PLT in µl`)/1000, 
    del5q = check_cols("del(5q)", "del(5)(q", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    del7 = check_cols("-7", "-7", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    plus8 = check_cols("+8", "+8", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    delY = check_cols("-Y", "-Y", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    del20q = check_cols("del(20q)", "del(20)(q", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    del7q = check_cols("del(7)(q", "del(7)(q", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    del11q = check_cols("del(11)(q", "del(11)(q", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    complex = check_cols("complex", "complex", `C Diagnosis`, `Karyotype 1`, `Karyotype 2`, `Karyotype 3`, `Karyotype 4`),
    MLL_PTD = as.numeric(NA),
    CYTO_IPSSR = case_when(
        complex == 1 ~ "Very Poor",
        del7 == 1 ~ "Poor",
        del7q == 1 | plus8 == 1 ~ "Intermediate",
        delY == 1 | del11q == 1 ~ "Very Good",
        TRUE ~ "Good"
    ))

hershberger_os <- os %>%
    mutate(ID = Exam_Array,
    OS_YEARS = OS_days/365.25,
    OS_STATUS = Censor_OS)

hershberger_mut <- mutations %>%
    mutate(across(-1, 
        ~ ifelse(. == "POSITIVE", 1,
        ifelse(. == "NEGATIVE", 0, NA)))) %>%
    left_join(tet2bi %>% select(ID, TET2bi), by = "ID") %>%
    mutate(TET2mono = ifelse(TET2 == 1 & TET2bi == 0, 1, 0))

hershberger_full <- hershberger %>%
    left_join(hershberger_os, by = "ID") %>%
    left_join(hershberger_mut, by = "ID") %>%
    mutate(TP53maxvaf = ifelse(WHO == "MDS-TP53", 0.7, 0),
    del17_17p = 0, 
    TP53loh = NA, 
    TET2mono = ifelse(TET2 == 1 & TET2bi == 0, 1, 0))

ipssm_process <- IPSSMprocess(hershberger_full)
ipssm_res <- IPSSMmain(ipssm_process)
ipssm_annot <- IPSSMannotate(ipssm_res)


classifySamples <- function(df){
    new_class <- case_when(
    df$complex == 1     ~ "Complex",
    df$del5q == 1      ~ "del5q-IB",
    df$SF3B1 == 1      ~ "SF3B1-IB", 
    df$EZH2 == 1        ~ "EZH2",
    df$TET2bi == 1      ~ "TET2-bi",
    df$del7 == 1       ~ "-7",
    df$STAG2 == 1       ~ "STAG2",
    df$BM_BLAST <= 5    ~ "MDS-LB",
    df$BM_BLAST > 10    ~ "MDS-IB2",
    df$BM_BLAST > 5 & df$BM_BLAST <= 10 ~ "MDS-IB1",
    TRUE                ~ "Other" 
  )
  factor(new_class, levels = c("EZH2",  "TET2-bi", "-7", "STAG2", "del5q-IB", "SF3B1-IB", "Complex",
      "MDS-LB", "MDS-IB1", "MDS-IB2"))
}


hershberger_full <- hershberger_full %>%
    left_join(ipssm_annot %>% 
    mutate(IPSSM = IPSSMcat_mean,
    IPSSM_SCORE = IPSSMscore_mean) %>%
    select(ID, IPSSM, IPSSM_SCORE), by = "ID") %>%
    mutate(sub_group = classifySamples(.))

hersh_all_mds <- filter(hershberger_full, grepl("MDS", WHO))
hersh_mds <- filter(hersh_all_mds, WHO %in% c("MDS-IB1", "MDS-IB2", "MDS-LB"))

save(hershberger_full, hersh_all_mds, file = "results/hershberger/hershberger_full.Rdata")
save(hersh_mds, file = "results/hershberger/hershberger_mds.Rdata")


## Create Summarized Experiment
counts <- read_table("data/hershberger/Counts.txt", skip = 1, col_names = FALSE)
header <- read_table("data/hershberger/Counts.txt", n_max = 1, col_names = FALSE)

counts_mat <- counts[, -1] %>% data.matrix()
colnames(counts_mat) <- header[1, ]
rownames(counts_mat) <- counts$X1


cData_hersh <- hershberger_full %>%
    filter(ID %in% colnames(counts_mat)) %>%
    column_to_rownames("ID")

## Add gene annotation
edb <- EnsDb.Hsapiens.v86
genes_coords <- genes(edb, filter = ~ symbol %in% rownames(counts_mat))
genes_coords_tb <- mcols(genes_coords) %>%
    as_tibble()
unique_coords <- group_by(genes_coords_tb, symbol)  %>%
    mutate(order = ifelse(gene_biotype == "protein_coding", 1, 2)) %>%
    arrange(symbol, order) %>%
    slice_head(n = 1)

genes_coords_sel <- genes_coords[match(unique_coords$gene_id, genes_coords$gene_id)]
names(genes_coords_sel) <- genes_coords_sel$symbol


final_ranges <- GRanges(
    seqnames = rep("*", nrow(counts_mat)),
    ranges = IRanges(start = 0, end = 0),
    gene_id = "",
    gene_name = rownames(counts_mat),
    gene_biotype = "",
    seq_coord_system = "chromosome", 
    symbol = rownames(counts_mat),
    entrezid = list("")
)
names(final_ranges) <- rownames(counts_mat)
final_ranges[names(genes_coords_sel)] <- genes_coords_sel


se_hersh <- SummarizedExperiment(assays = list(counts = counts_mat), 
    colData = cData_hersh[colnames(counts_mat), ],
    rowRanges = final_ranges)
save(se_hersh, file = "results/hershberger/rnaseq_hersh.Rdata")