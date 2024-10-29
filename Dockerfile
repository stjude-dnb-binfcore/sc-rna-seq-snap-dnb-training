FROM rocker/rstudio:4.4.0

# Set environment variables
ENV PATH=/usr/lib/rstudio-server/bin:${PATH}
ENV R_LIBS_USER=${HOME}/R/my-rstudio/4.4.0
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Install base dependencies and locales
RUN apt-get update 
RUN apt-get install -y --no-install-recommends \
        locales \
        aptitude \
        gnupg \
        software-properties-common \
        dirmngr 
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen 
RUN locale-gen en_US.utf8 
RUN /usr/sbin/update-locale LANG=en_US.UTF-8 
RUN apt-get clean 
RUN rm -rf /var/lib/apt/lists/*

# Install R and dependencies
RUN apt-get update && \
    apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libcairo2-dev \
        libxt-dev \
        libpng-dev \
        libfreetype6-dev \
        libtiff5-dev \
        libjpeg-dev \
        libgsl-dev \
        librsvg2-dev \
        libnode-dev \
        libv8-dev \
        software-properties-common \
        libharfbuzz-dev \
        librsvg2-dev \
        libproj-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
        libfribidi-dev \
        python3-pip \
        python3-dev \
        cmake \
        bedtools \
        openjdk-8-jre-headless \
        g++ \
        libopenblas-base \
        liblapack3 \
        libgeos-dev \
        pkg-config \
        jags \
        git-all \
        libudunits2-dev \
        libmagick++-dev \
        libgdal-dev \
        libhdf5-dev \
        libpoppler-cpp-dev \
        libfftw3-dev \
        libglpk-dev \
        libgmp-dev \
        libgit2-dev \
        curl \
        build-essential \
        autoconf \
        automake \
        flex \
        bison && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

# Install Python packages
RUN pip3 install MACS2 leidenalg

# Enable RStudio Copilot
RUN echo "copilot-enabled=1" >> /etc/rstudio/rsession.conf

# Install R packages
RUN Rscript -e 'install.packages("usethis")' && \
    Rscript -e 'install.packages(c("BiocManager", "devtools", "remotes"))' 

RUN Rscript -e 'BiocManager::install(c("miQC", "scater", "scDblFinder", "SingleCellExperiment", "celldex", "SingleR", "infercnv"))' 

RUN Rscript -e 'devtools::install_github("igordot/scooter")' 
RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc")' 
RUN Rscript -e 'devtools::install_github("kharchenkolab/numbat")' 

RUN Rscript -e 'install.packages(c("clustree", "cowplot", "data.table"))'
RUN Rscript -e 'install.packages(c("flexmix", "flextable", "forcats"))'
RUN Rscript -e 'install.packages(c("fs", "future", "GGally"))'
RUN Rscript -e 'install.packages(c("ggh4x", "ggplot2", "ggpmisc"))'
RUN Rscript -e 'install.packages(c("ggrepel", "ggthemes", "grid"))'
RUN Rscript -e 'install.packages(c("harmony", "igraph", "irlba"))'
RUN Rscript -e 'install.packages(c("knitr", "leiden", "optparse"))'
RUN Rscript -e 'install.packages(c("patchwork", "purrr", "RColorBrewer"))'
RUN Rscript -e 'install.packages(c("remotes", "reshape2", "rliger"))'
RUN Rscript -e 'install.packages(c("rlist", "R.utils", "SeuratObject"))'
RUN Rscript -e 'install.packages(c("shiny", "SoupX", "stringr"))'
RUN Rscript -e 'install.packages(c("tidytext", "tidyverse", "tinytex", "yaml"))' 

# Seurat
RUN Rscript -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'

# Set the entrypoint to rserver
ENTRYPOINT ["/init"]
