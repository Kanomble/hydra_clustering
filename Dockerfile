# DOCKERFILE for constructing the deseq2 container - base image of kanomble/deseq2
FROM rocker/rstudio:4.1.2
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk40 \
    libgit2-28 \
  && apt-get clean all && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN Rscript -e "install.packages(c('devtools','rmarkdown', 'tidyverse'));"
RUN Rscript -e "install.packages(c('BiocManager'));"
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("DEGreport")'
RUN Rscript -e "install.packages(c('pheatmap','RColorBrewer','PoiClaClu','ggbeeswarm'));"
RUN Rscript -e "install.packages(c('shiny'));"
RUN Rscript -e "install.packages(c('reshape2','viridis','scales','ggdendro','gridExtra'));"