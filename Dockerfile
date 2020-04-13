FROM docker.io/ubuntu:18.04
FROM openjdk:11
FROM rocker/rstudio
FROM rocker/verse
MAINTAINER hajk-georg.drost@tuebingen.mpg.de
RUN ["apt-get", "update"]
RUN ["apt-get", "-y", "install", "apt-utils"]
RUN ["apt-get", "-y", "install", "gcc"]
RUN ["apt-get", "-y", "install", "python3"]
RUN ["apt-get", "-y", "install", "perl"]
RUN ["apt-get", "-y", "install", "make"]
RUN ["apt-get", "-y", "install", "sudo"]
RUN ["apt-get", "-y", "install", "wget"]
RUN ["apt-get", "-y", "install", "genometools"]
RUN ["apt-get", "-y", "install", "git"]
RUN ["apt-get", "-y", "install", "autoconf"]
RUN ["apt-get", "-y", "install", "g++"]
RUN ["apt-get", "-y", "install", "ncbi-blast+"]
RUN ["apt-get", "-y", "install", "build-essential"]
RUN ["apt-get", "-y", "install", "libcurl4-gnutls-dev"]
RUN ["apt-get", "-y", "install", "libxml2-dev"]
RUN ["apt-get", "-y", "install", "libssl-dev"]
RUN ["apt-get", "-y", "install", "libpq-dev"]
RUN ["apt-get", "-y", "install", "software-properties-common"]
RUN sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN ["apt-get", "-y", "install", "r-base"]
RUN apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
WORKDIR "/app"
RUN mkdir software_downloads \
  && cd software_downloads \
  && wget http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz \
  && tar xf hmmer-3.3.tar.gz \
  && cd hmmer-3.3 \
  && ./configure \
  && make \
  && make check \
  && sudo make install \
  && cd ..
RUN cd software_downloads \
  && wget https://github.com/torognes/vsearch/archive/v2.14.2.tar.gz \
  && tar xzf v2.14.2.tar.gz \
  && cd vsearch-2.14.2 \
  && sudo ./autogen.sh \
  && sudo ./configure \
  && sudo make \
  && sudo make install \
  && cd ..
RUN cd software_downloads \
  && wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan.pl.gz \
  && gunzip dfamscan.pl.gz \
  && sudo cp dfamscan.pl /usr/local/bin/dfamscan.pl \
  && cd ..
RUN cd software_downloads \
  && wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz \
  && gunzip usearch11.0.667_i86linux32.gz \
  && chmod +x usearch11.0.667_i86linux32 \
  && sudo mv usearch11.0.667_i86linux32 usearch \
  && sudo cp usearch /usr/local/bin/usearch \
  && cd ..
RUN sudo R -e "install.packages('devtools')"
RUN sudo R -e "install.packages('tidyverse')"
RUN sudo R -e "install.packages('BiocManager')"
RUN sudo R -e "BiocManager::install()"
RUN sudo R -e "BiocManager::install(c('rtracklayer', 'GenomicFeatures', 'GenomicRanges', 'GenomeInfoDb', 'biomaRt', 'Biostrings'))"
RUN sudo R -e "install.packages(c('tidyverse', 'data.table', 'seqinr', 'biomartr', 'ape', 'dtplyr', 'devtools'))"
RUN sudo R -e "devtools::install_github('HajkD/metablastr', build_vignettes = TRUE, dependencies = TRUE)"
RUN sudo R -e "install.packages(c('BSDA', 'ggrepel', 'gridExtra'))"
RUN wget https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz
RUN sudo R -e "install.packages('latticeExtra_0.6-28.tar.gz', type = 'source')"
RUN sudo R -e "install.packages('survival')"
RUN sudo R -e "BiocManager::install('ggbio')"
RUN sudo R -e "devtools::install_github('HajkD/LTRpred')"
LABEL version="1.1.0"
VOLUME ["/data"]
EXPOSE 80
CMD ["/bin/bash"]