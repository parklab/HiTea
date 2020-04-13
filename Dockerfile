# Source
FROM ubuntu:19.10

# Maintainer
# MAINTAINER Dhawal Jain (dhawal.sjain@gmail.com)

# general updates & installing necessary Linux components
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update -y 
RUN apt-get install -y \
    bzip2 \
    build-essential \
    gcc \
    gcc-multilib \
    apt-utils \
    expat \
    libexpat-dev \ 
    make \
    man \
    parallel \
    perl \
    unzip \
    libcurl4-openssl-dev \
    libssl-dev \
    r-base \
    r-base-dev \
    time \
    vim \
    wget \
    bedtools

# installing conda
RUN wget https://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh && bash Anaconda3-4.4.0-Linux-x86_64.sh -p /anaconda3 -b
ENV PATH=/anaconda3/bin:$PATH
RUN conda update -y conda \
    && rm Anaconda3-4.4.0-Linux-x86_64.sh


# installing R libraries from CRAN
RUN R -e "install.packages(c('data.table', 'BiocManager','rmarkdown', 'knitr', 'ggplot2','circlize','kableExtra','DT'), \
          repos='http://cran.us.r-project.org/' ,clean = TRUE )"

# installing R libraries from Bioconductor
RUN R -e "BiocManager::install(pkgs = c('GenomicRanges','EnrichedHeatmap'),ask=FALSE,clean = TRUE)"


# get pandoc for rmarkdown report generation (optional)
RUN apt install -y texlive-latex-base pandoc texlive-fonts-recommended


# set work directory
WORKDIR /usr/local/bin


# bwa 0.7.17
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz \
        && tar -xzf v0.7.17.tar.gz \
        && cd bwa-0.7.17 \
        && make \
        && cd .. \
        && ln -s bwa-0.7.17 bwa


# samtools 1.7
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
         && tar -xjf samtools-1.7.tar.bz2 \
         && cd samtools-1.7 \
         && make \
         && cd .. \
         && ln -s samtools-1.7 samtools 

# set path
ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH

# pairtools 0.3.0
RUN conda config --add channels conda-forge \
          && conda config --add channels bioconda \
          && conda install -c conda-forge -c bioconda pairtools

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8


# cleanup
RUN apt-get autoclean && apt-get clean
RUN rm samtools-1.7.tar.bz2 && rm v0.7.17.tar.gz

# hitea 0.1.01
RUN wget https://github.com/parklab/HiTea/archive/0.1.01.tar.gz \
        && tar -xzf 0.1.01.tar.gz \
        && rm 0.1.01.tar.gz    


# Copy script and data directory
RUN mv HiTea*/* .
RUN rm -d HiTea*
COPY test.sh .

#COPY hitea .
#ADD examples /usr/local/bin/examples
#ADD hg19 /usr/local/bin/hg19
#ADD hg38 /usr/local/bin/hg38
#ADD src /usr/local/bin/src

RUN chmod +x hitea


# default command
CMD ["hitea"]