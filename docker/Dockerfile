############################################################
# Part 1: Base image

FROM ubuntu:22.04
LABEL maintainer="YangHang <hangyang@genetics.ac.cn>"

############################################################
# Part 2: System apps

WORKDIR /tmp

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        git wget curl unzip vim emacs sudo \
        build-essential software-properties-common gdebi-core llvm \
        pandoc pandoc-citeproc \
        apache2 apache2-utils \
        perl ruby ruby-dev nodejs npm \
        tzdata zlib1g-dev libcurl4-gnutls-dev libcairo2-dev libxt-dev \
        libssl-dev libssh2-1-dev libmpfr-dev libxml2-dev libharfbuzz-dev \
        libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libbz2-dev libreadline-dev libsqlite3-dev \
        libncurses5-dev libncursesw5-dev xz-utils tk-dev \
        libffi-dev liblzma-dev

# Python
# default: system Python, also install 2.7 in pyenv for some apps
RUN apt-get update && \
    apt-get install -y python3 python3-pip python-is-python3 \
            python3-chardet python3-pygraphviz && \
    pip install --upgrade pip && \
    pip install jinja2 requests chardet "dask<=2023.9.1" "distributed==2023.9.1"

ENV PYENV_ROOT /root/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH
RUN git clone https://github.com/pyenv/pyenv.git $PYENV_ROOT && \
    pyenv install 2.7 && \
    pyenv global  2.7 && \
    pip install --upgrade pip && \
    pip install numpy && \
    pyenv global system

# R
RUN apt-get install -y r-base && \
    R -e "install.packages(c('shiny', 'shinyWidgets', 'ggplot2', 'scales', 'reshape2', 'stringr', 'tidyr', 'grid', 'gridExtra', 'cowplot', 'magrittr', 'RColorBrewer', 'RSQLite', 'shinycssloaders', 'shinythemes', 'SeuratObject', 'BiocManager', 'xml2', 'ragg', 'devtools', 'patchwork'), repos='https://cloud.r-project.org/')" && \
    wget -q https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb && \
    gdebi -n shiny-server-1.5.20.1002-amd64.deb

# Vue
RUN npm install -g @vue/cli

############################################################
# Part 3: Bioinfo apps

RUN apt-get install -y \
        fastqc trim-galore hisat2 bowtie2 stringtie samtools \
        picard-tools python3-deeptools bedtools gffread ncbi-blast+
RUN R -e "BiocManager::install(c('S4Vectors', 'GenomeInfoDb'))"

# UCSC tools
ARG TOOLS="bedClip bedGraphToBigWig bigWigAverageOverBed genePredToBed genePredToGtf bedToGenePred gtfToGenePred gff3ToGenePred"
RUN mkdir -p /usr/local/bin && \
    (cd /usr/local/bin && \
    for tool in $TOOLS; do \
        wget -q "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$tool" && \
        chmod +x "$tool"; \
    done)

# sra-toolkit
RUN wget -q --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -C /usr/local/share/ -xzf sratoolkit.tar.gz && \
    ln -s /usr/local/share/sratoolkit.3.0.7-ubuntu64/bin/* /usr/local/bin/

# MACS2
RUN pip install rseqc MACS2

# bismark
RUN wget -q --output-document Bismark.zip https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.24.1.zip && \
    unzip Bismark.zip -d /usr/local/share/ && \
    find /usr/local/share/Bismark-0.24.1/ -type f -executable -exec ln -s {} /usr/local/bin/ \;

# ROSE
RUN wget -q --output-document ROSE-master.zip https://github.com/stjude/ROSE/archive/refs/heads/master.zip && \
    unzip ROSE-master.zip -d /usr/local/share/ && \
    ln -s /usr/local/share/ROSE-master/lib/* /usr/local/bin/ && \
    ln -s /usr/local/share/ROSE-master/bin/* /usr/local/bin/

# ANANSE
RUN pip install git+https://github.com/vanheeringen-lab/ANANSE

# SequenceServer
RUN gem install sequenceserver
# RUN mkdir -p ~/.sequenceserver && \
#     touch ~/.sequenceserver/asked_to_join && \
#     mkdir -p /db/blastdb && \
#     cd /db/blastdb && \
#     echo '>test \nATGCATGC \n' > test.fa && \
#     makeblastdb -in test.fa -dbtype nucl -parse_seqids

# JBrowse
RUN wget -q https://github.com/GMOD/jbrowse/releases/download/1.16.11-release/JBrowse-1.16.11.zip && \
    unzip JBrowse-1.16.11.zip -d /var/www/html/ && \
    mv /var/www/html/JBrowse-1.16.11 /var/www/html/jbrowse_default && \
    /var/www/html/jbrowse_default/setup.sh

# prepare miniODP data links and working dir
RUN ln -s /mnt/miniENCODE/miniODP/miniodp /var/www/html/miniodp && \
    ln -s /mnt/miniENCODE/miniODP/jbrowse /var/www/html/jbrowse && \
    ln -s /mnt/miniENCODE/miniODP/mODP_zebrafish /srv/shiny-server/mODP_zebrafish && \
    mkdir -p /db/blastdb && ln -s /mnt/miniENCODE/miniODP/blastdb /db/blastdb

WORKDIR /mnt

# Expose necessary ports
EXPOSE 80 3838 4040

# Clean up
RUN rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* && apt-get clean
