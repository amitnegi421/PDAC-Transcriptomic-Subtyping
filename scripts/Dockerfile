FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN R -e "BiocManager::install(c('TCGAbiolinks', 'SummarizedExperiment', 'DESeq2', 'clusterProfiler', 'GEOquery'))" && \
    R -e "install.packages(c('tidyverse', 'survival', 'survminer', 'pheatmap', 'msigdbr', 'tidyestimate', 'enrichR'))"

WORKDIR /workspace
COPY scripts/ /workspace/scripts/