FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get install -y build-essential unzip wget 
RUN apt-get install -y --no-install-recommends openjdk-11-jre
RUN apt-get install -y time locales 
RUN apt-get install -y make python3-pip \
  python3-dev \
  cpanminus 
RUN apt-get install -y  bwa



WORKDIR /usr/local/

RUN wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && \
  unzip fastqc_v0.11.5.zip && \
  rm fastqc_v0.11.5.zip && \
  cd FastQC && \
  chmod 755 fastqc && \
  ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

RUN pip install multiqc

ENV LC_ALL C
ENV PATH=/usr/local/bin:$PATH

RUN apt-get install -y libncurses5-dev libbz2-dev liblzma-dev
# samtools
ARG SAMTOOLSVER=1.16.1
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install 

RUN apt-get install -y python

# GATK
RUN wget 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip' && \
  unzip gatk-4.3.0.0.zip

ENV PATH=/usr/local/gatk-4.3.0.0/:${PATH}


#Freebayes

RUN apt-get install -y bc 
RUN apt-get install -y parallel meson 
RUN apt-get install -y ninja-build 
RUN apt-get install -y libvcflib-tools vcftools

RUN wget 'https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-src.tar.gz' && \
  tar -zxvf freebayes-1.3.6-src.tar.gz && \
  cd freebayes && \
  meson build/ --buildtype release && \
  cd build && ninja 

ENV PATH=/usr/local/freebayes/build/:/usr/local/freebayes/scripts/:${PATH}


# Cleanup apt package lists to save space
RUN rm -rf /var/lib/apt/lists/*

