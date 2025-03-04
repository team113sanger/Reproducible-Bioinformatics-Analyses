FROM rocker/tidyverse:4.2.0

# Install additional packages from CRAN
RUN install2.r --error \
    rmarkdown \
    shiny \
    knitr

# Install packages from Bioconductor
RUN R -e "BiocManager::install(c('DESeq2', 'edgeR'))"
