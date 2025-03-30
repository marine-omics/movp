FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential unzip wget openjdk-17-jre time locales \
  make python-is-python3 python3-pip \
  python3-dev \
  cpanminus bwa \
  libncurses5-dev libbz2-dev liblzma-dev \
  bc parallel meson ninja-build libvcflib-tools vcftools zlib1g pkg-config cmake


WORKDIR /usr/local/

RUN wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && \
  unzip fastqc_v0.11.5.zip && \
  rm fastqc_v0.11.5.zip && \
  cd FastQC && \
  chmod 755 fastqc && \
  ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

RUN pip install multiqc

ENV LC_ALL=C
ENV PATH=/usr/local/bin:$PATH


# samtools
ARG SAMTOOLSVER=1.21
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install && rm -rf /usr/local/samtools-${SAMTOOLSVER}

#bcftools
ARG BCFTOOLSVER=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLSVER}/bcftools-${BCFTOOLSVER}.tar.bz2 && \
  tar -xjf bcftools-${BCFTOOLSVER}.tar.bz2 && \
 rm bcftools-${BCFTOOLSVER}.tar.bz2 && \
 cd bcftools-${BCFTOOLSVER} && \
 ./configure && \
 make && \
 make install && rm -rf /usr/local/bcftools-${BCFTOOLSVER}



# GATK
RUN wget 'https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip' && \
  unzip gatk-4.6.1.0.zip && rm gatk-4.6.1.0.zip

ENV PATH=/usr/local/gatk-4.6.1.0/:${PATH}


RUN wget 'https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-src.tar.gz' && \
  tar -zxvf freebayes-1.3.6-src.tar.gz && \
  cd freebayes && \
  meson build/ --buildtype release && \
  cd build && ninja && rm /usr/local/freebayes-1.3.6-src.tar.gz

ENV PATH=/usr/local/freebayes/build/:/usr/local/freebayes/scripts/:${PATH}


# Fastp
WORKDIR /usr/local/bin/
RUN wget http://opengene.org/fastp/fastp.0.23.1 && \
    mv fastp.0.23.1 fastp &&\
    chmod a+x ./fastp

#htslib
ARG HTSLIBVER=1.21
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIBVER}/htslib-${HTSLIBVER}.tar.bz2 && \
  tar -xjf htslib-${HTSLIBVER}.tar.bz2 && \
 rm htslib-${HTSLIBVER}.tar.bz2 && \
 cd htslib-${HTSLIBVER} && \
 ./configure && \
 make && \
 make install && rm -rf /usr/local/htslib-${HTSLIBVER}

#
# ------ Installing vcflib / vcfwave.  Installs but vcfwave seems buggy as of 24Feb 2023 ----- #
#
# RUN apt-get update && apt-get install -y --no-install-recommends \
#   git cmake libhts-dev libtabixpp-dev libtabixpp0 \
#   pybind11-dev libcurl4-openssl-dev python-dev \
#   libomp-dev

# WORKDIR /usr/local/

# RUN git clone --recursive https://github.com/vcflib/vcflib.git && \
#   cd vcflib && \
#   mkdir -p build && cd build && \
#   cmake -DCMAKE_BUILD_TYPE=Release -DZIG=OFF .. && \
#   cmake --build . && \
#   cmake --install .

#ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/vcflib/build/contrib/WFA2-lib/

# Cleanup apt package lists to save space
RUN rm -rf /var/lib/apt/lists/*

