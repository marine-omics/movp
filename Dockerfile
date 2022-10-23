FROM ubuntu:22.04


RUN apt-get update \
	&& apt-get install -y build-essential unzip wget default-jre \
	python3-pip

# FASTQC
ENV DST=/tmp
ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v0.11.5.zip

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  unzip - $DST/$ZIP -d $DST && \
  rm $DST/$ZIP && \
  cd $DST/FastQC && \
  chmod 755 fastqc && \
  ln -s $DST/FastQC/fastqc /usr/local/bin/fastqc

RUN pip install multiqc

ENV PATH /usr/local/bin:$PATH

