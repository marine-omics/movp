FROM ubuntu:20.04


RUN apt-get update
RUN apt-get install -y build-essential unzip wget 
RUN apt-get install -y --no-install-recommends openjdk-11-jre
RUN apt-get install -y time locales 
RUN apt-get install -y make python3-pip \
  python3-dev \
  cpanminus \
  && rm -rf /var/lib/apt/lists/*


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


