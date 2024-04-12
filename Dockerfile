FROM rocker/shiny-verse:4.3.3

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
    libbz2-dev \
    liblzma-dev \
    g++ \
    gcc \
    libglpk-dev \
    software-properties-common && add-apt-repository -y ppa:git-core/ppa

# main targets
# Seurat => V5 at this point
RUN R -e "install.packages(c('plotly', \
    'markdown', \
    'shinyjs', \
    'ggpubr', \
    'shinyhelper', \
    'shinycssloaders', \
    'dplyr', \
    'bslib', \
    'Seurat'))"

RUN R -e "remotes::install_github('bnprks/BPCells')"

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