# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install basic R packages
RUN R -e 'BiocManager::install(c("cowplot", "tidyverse", "parallel", "caret", "e1071", "GGally"))'

# Install clustering R packages
RUN R -e 'BiocManager::install(c("FactoMineR", "missMDA","FactoInvestigate", "clustMixType", "mclust"))'

# Install survival R packages
RUN R -e 'BiocManager::install(c("survival", "survminer", "cmprsk"))'

# Install graphical packages (v.1.1)
RUN R -e 'BiocManager::install(c("pheatmap"))'

# Install tsne (v.1.2)
RUN R -e 'BiocManager::install(c("Rtsne"))'

# Install ipssm (v.1.3)
RUN R -e 'devtools::install_github("papaemmelab/ipssm",ref = "main")'

# Install umap (v.1.4)
RUN R -e 'BiocManager::install(c("umap"))'

# Install GO related packages (v.1.5)
RUN R -e 'BiocManager::install(c("GOfuncR", "Homo.sapiens"))'

# Install ca for correspondency analysis (v.1.6)
RUN R -e 'BiocManager::install(c("ca"))'

# Install forestploter and RNAseq  for forest plots (v.1.7)
RUN R -e 'BiocManager::install(c("forestploter", "DESeq2", "edgeR", "limma", "NetActivity", "GSVA"))'

# Install genes annotation (v.1.8)
RUN R -e 'BiocManager::install(c("EnsDb.Hsapiens.v86"))'