FROM ubuntu:18.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 


RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python \
    python-pip

RUN wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    tar jxvf htslib-1.7.tar.bz2 && \
    cd htslib-1.7 && \
    ./configure && \
    make && \
    make install

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN pip install annot_utils==0.3.0
RUN pip install pysam==0.15
RUN pip install numpy==1.15.1
RUN pip install scipy==1.1.0
RUN pip install statistics==1.0.3.5

COPY EBbreak_modified.tar.gz /

RUN tar zxvf EBbreak_modified.tar.gz && cd EBbreak_modified &&  pip install . 

