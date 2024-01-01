# Part 1 : Base image
FROM openanalytics/r-base

LABEL maintainer "YangHang <hangyang@genetics.ac.cn>"

# Part 2 : system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    git \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libmpfr-dev \
    build-essential \
    libxml2-dev \
    libgeos-dev \
    gfortran \
    libgsl-dev

# Part 3 : The R package will be used
RUN R -e "install.packages(c('shiny', 'shinyWidgets', 'ggplot2', 'scales', 'reshape2', 'stringr', 'tidyr', 'grid', 'gridExtra', 'cowplot', 'magrittr', 'RColorBrewer', 'RSQLite', 'dplyr', 'Seurat', 'shinycssloaders', 'devtools', 'BiocManager', 'uwot', 'RcppParallel', 'shinythemes', 'visNetwork'), repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('S4Vectors', 'chromVAR'), library = '/usr/local/lib/R/site-library')"
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories(), library = '/usr/local/lib/R/site-library')"
RUN R -e "remotes::install_github('tanaylab/shaman')"

# Part 4 : Copy the app to the image
RUN mkdir -p /mnt/mODPbase/

# Part 5 : Copy Rprofile.site
RUN echo -e 'local({\n options(shiny.port = 3838, shiny.host = "0.0.0.0")\n})' > /usr/lib/R/etc/Rprofile.site

# Part 6 : Expose port, must equal to the port in Rprofile.site
EXPOSE 3838

# Part 7 : Command to run the shiny app
CMD ["R", "-e", "shiny::runApp('/mnt/mODPbase')"]
