# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install basic R packages
RUN R -e 'BiocManager::install(c("cowplot", "tidyverse", "parallel", "caret", "e1071", "GGally"))'

# Install clustering R packages
RUN R -e 'BiocManager::install(c("FactoMineR", "missMDA","FactoInvestigate", "clustMixType", "mclust"))'

# Install survival R packages
RUN R -e 'BiocManager::install(c("survival", "survminer", "cmprsk"))'
