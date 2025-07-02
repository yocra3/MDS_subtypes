#' ---------------------------
#'
#' Purpose of script:
#'
#' Define GO terms to use for predicting risk
#'
#' ---------------------------
#'
#' Notes:
#'
#' Docker command:   
#' docker run -it -v $PWD:$PWD -w $PWD mds_subtypes_rsession:1.5 R
#'
#' ---------------------------

## Load libraries 
library(GOfuncR)
library(tidyverse)
library(parallel)
library(caret)
library(survival)

## Load new GO graph
godir <-  "data/GO_terms/"
term <- read.table(paste0(godir, "/term.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)
graph <- read.table(paste0(godir, "/graph_path.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)

## Get genes present in data
maf <- read_tsv("./data/IPSSMol/maf.tsv")
genes <- unique(maf$GENE)

## Get all GO terms
all_gos <- get_child_nodes("GO:0008150", term, graph)

## Remove obsolete terms
term_filt <- subset(term, V5 == 0)
all_gos <- subset(all_gos, child_go_id %in% term_filt$V4)
genes_pairs <- get_anno_genes(go_ids = all_gos$child_go_id, term_df = term, graph_path_df = graph)

tab <- left_join(mutate(all_gos, go_id = child_go_id), genes_pairs, by = "go_id") %>%
  as_tibble() %>%
  filter(!is.na(gene))

tab$PathwayID <- tab$go_id

## Compute proportion of genes present in GO
gos_gene_tab <- tab %>%
  group_by(go_id) %>%
  summarize(prop = mean(gene %in% genes),
            n_total = n(),
            n_genes = sum(gene %in% genes))

## Select GO terms with at least 3 genes present and minimum 3% of overlap
gos_gene_sel <- subset(gos_gene_tab, n_genes >= 3 & prop > 0.03)
pres_genes <- sum(genes %in% subset(tab, go_id %in% gos_gene_sel$go_id)$gene) ## All genes present in this subset

gene_tab <- table(maf$GENE)
gene_tab[!names(gene_tab) %in% subset(tab, go_id %in% gos_gene_sel$go_id)$gene]


## Select GO terms in order of maximum overlap with selection. 
sel_genes <- c()
sel_gos <- c()
i <- 1
gos_gene_ord <- arrange(gos_gene_sel, desc(prop))
while(length(sel_genes) < pres_genes){
  message("Iteration: ", i, " N genes: ", length(sel_genes))
  sel_go <- gos_gene_ord[i, ]$go_id
  new_genes <- genes[genes %in% subset(tab, go_id %in% sel_go)$gene]
  if (!all(new_genes %in% sel_genes)){
    sel_gos <- c(sel_gos, sel_go)
    sel_genes <- union(sel_genes, new_genes)
  }
  i <- i + 1
}
tab_final <- filter(tab, go_id %in% sel_gos & gene %in% sel_genes)
mis_genes <- genes[which(!genes %in% tab_final$gene)]
## 6 discarded with this approach.

length(unique(tab_final$PathwayID))
# [1] 66

gos_list <- split(tab_final$gene, tab_final$PathwayID)

gene_mat <- matrix(0, length(gos_list), length(unique( tab_final$gene)),
 dimnames = list(names(gos_list), unique(tab_final$gene)))
for (go in names(gos_list)) 
    gene_mat[go, gos_list[[go]]] <- 1

gene_d <- dist(gene_mat, "binary")
gene_dmat <- as.matrix(gene_d)
gene_sim <- gene_dmat 
diag(gene_sim) <- 1

## Select GO terms sharing < 50% of genes 
## From these groups, retain the largest GO 
sel_paths <- names(gos_list)

gene_dloop <- gene_d
gene_dmat_loop <- gene_dmat
path_cls <- cutree(hclust(gene_dloop), h = 0.5)
round <- 1
while(length(sel_paths) != length(unique(path_cls))){
    message("Round ", round)
  sel_paths <- sapply(unique(path_cls), function(cl){
    paths <- rownames( gene_dmat_loop)[path_cls == cl]
    cl_gos <- gos_list[paths]
    names(cl_gos)[which.max(lengths(cl_gos))]
  })
  gene_dmat_loop <- gene_dmat_loop[sel_paths, sel_paths]
  gene_dloop <- as.dist(gene_dmat_loop)
  path_cls <- cutree(hclust(gene_dloop), h = 0.5)
  round <- round + 1
}

tab_out <- subset(tab_final, PathwayID %in% sel_paths)
length(unique(tab_out$go_id))
## 59 pathways
length(unique(tab_out$gene))
## 115 genes

save(tab_final, tab_out, file = "results/mutations/go_gene_map.Rdata")

write.table(tab_final[, c("PathwayID", "gene")], file = "results/mutations/go_gene_map_unfiltered.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


write.table(tab_out[, c("PathwayID", "gene")], file = "results/mutations/go_gene_map_filtered.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
