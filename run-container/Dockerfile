FROM rocker/rstudio:4.4.0

# Set environment variables
#ENV PATH=/usr/lib/rstudio-server/bin:${PATH}
ENV PATH=/opt/cellranger-8.0.1/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/cellranger-8.0.1/lib:${PATH}
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
        python3.9 \
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
        bison \
        texlive \
        pandoc \
        fastqc \
        gdebi-core && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install rsession
RUN wget -P /tmp https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.09.1-394-amd64.deb && \
    gdebi -n /tmp/rstudio-server-2024.09.1-394-amd64.deb

# Install CellRanger v8.0.1
RUN curl -o /tmp/cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1731127395&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=V6lD7p-GGp2O~XSxRC3K9olOPKUYtLGB0oUH845lEIzFgYdBsNoMVhEuNm2u1p~ndCBxGLYM1~PCDLXt9aZM5qQUtWRbwfMqTmuCJ30JMa9UGZY99SWzeE7B3h-CD4Vk6PoR0xBwJPLH5uG4B-kLgNlXSAzSO~lQMjjupoP2tYD37~3YAf91RJilOmaa0t-PEMHMeTTvlgThl6M9zQ08Y2VuJquA7H-zDkblwUSmqexR8J22YF-CtvyYVm~KgZqHezi0r96r5~oeSICmsniBF5KBEI6I1ERBh5MNs2w-kx03V3fkGAOXwcbelfM3mBWtPMElMvekQuHj5bBqz1GvUA__"
RUN tar -xvzf /tmp/cellranger-8.0.1.tar.gz -C /opt

# Install pip for Python 3.9
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3

# Set default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

# Install Python packages
RUN pip3 install MACS2 leidenalg multiqc

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
RUN Rscript -e 'install.packages(c("knitr", "optparse"))'
RUN Rscript -e 'install.packages(c("patchwork", "purrr", "RColorBrewer"))'
RUN Rscript -e 'install.packages(c("remotes", "reshape2", "rliger"))'
RUN Rscript -e 'install.packages(c("rlist", "R.utils", "SeuratObject"))'
RUN Rscript -e 'install.packages(c("shiny", "SoupX", "stringr"))'
RUN Rscript -e 'install.packages(c("tidytext", "tidyverse", "tinytex", "yaml"))' 

RUN Rscript -e 'devtools::install_github("welch-lab/liger")'

# Seurat
RUN Rscript -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'

# Set the entrypoint to rserver
ENTRYPOINT ["/init"]
