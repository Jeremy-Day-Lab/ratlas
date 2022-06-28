FROM rocker/shiny-verse:4.1.0

LABEL maintainer="Lara Ianov <lianov@uab.edu>"
LABEL description="Dockerized Ratlas app dependencies from the Day lab - UAB"

#ignores conda install prompt from reticulate to keep log cleaner
COPY ./.Renviron /home/shiny/.Renviron

#custom server config:
COPY ./shiny_customized.config /etc/shiny-server/shiny-server.conf

RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    llvm-10 \
    libbz2-dev \
    liblzma-dev \
    g++ \
    gcc \
    libglpk-dev \
    software-properties-common && add-apt-repository -y ppa:git-core/ppa

# BioC packages

RUN R --no-restore --no-save -e "BiocManager::install(c('S4Vectors', 'SummarizedExperiment', \
'SingleCellExperiment', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', \
'Biobase', 'GenomeInfoDbData'))"

# additional suggestions (CRAN) + some of our previouslly needed packages

RUN R --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'Rfast2','hdf5r', 'remotes'))"

# main targets
RUN R -e "install.packages('plotly')"
RUN R -e "install.packages('markdown')"
RUN R -e "install.packages('shinyjs')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('Seurat')"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN R -e "install.packages('shinyhelper')"
RUN R -e "install.packages('shinycssloaders')"
RUN R -e "install.packages('dplyr')"
RUN R -e "setRepositories(ind=1:2); install.packages('Signac',  repos = c('https://cloud.r-project.org', 'https://bioconductor.org/packages/3.13/bioc'),dependencies = TRUE)"

# port container listens to (default R Shiny port)
EXPOSE 3838

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server
RUN sudo chown -R shiny:shiny /var/log/shiny-server/

# remove sample index and apps:
RUN sudo rm /srv/shiny-server/index.html
RUN sudo rm -rf /srv/shiny-server/sample-apps
#exec form for wildcards...
RUN ["/bin/bash", "-c", "sudo rm -rf /srv/shiny-server/{0*,1*}_*"]

# run app
CMD ["/usr/bin/shiny-server"]