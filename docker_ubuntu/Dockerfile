FROM ubuntu:20.04

LABEL maintainer "epil02@uw.edu"

ARG BIOPYTHON_VERSION=1.78
ARG MUMMER_VER="4.0.0"
ARG GAP2SEQ_VER="3.1"
ARG SAMTOOLS_VER="1.15"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=$PATH:/usr/local/bin

# Update package list and install necessary packages
RUN apt-get update && apt-get install -y --no-install-recommends \
  wget \
  git \
  bc \
  libncurses5-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  zlib1g-dev \
  libssl-dev \
  gcc \
  make \
  cmake \
  perl \
  python3 \
  python-is-python3 \
  python3-pip \
  bzip2 \
  gnuplot \
  ca-certificates \
  gawk \
  curl \
  sed \
  gnuplot \
  build-essential \
  unzip && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

# make path folder
RUN mkdir -p /usr/local/bin

# Install Biopython
RUN python3 -m pip install biopython==${BIOPYTHON_VERSION}

# download mummer
WORKDIR /tmp

RUN wget https://github.com/mummer4/mummer/releases/download/v${MUMMER_VER}rc1/mummer-${MUMMER_VER}rc1.tar.gz && \
  tar -xvf mummer-${MUMMER_VER}rc1.tar.gz && \
  rm mummer-${MUMMER_VER}rc1.tar.gz && \
  cd mummer-${MUMMER_VER}rc1 && \
  ./configure && \
  make && \
  make install && \
  ldconfig && \
  cd ..

RUN nucmer -h 

# download gap2seq
RUN wget https://github.com/rikuu/Gap2Seq/releases/download/v${GAP2SEQ_VER}/gap2seq-${GAP2SEQ_VER}-x64-linux.tar.bz2 && \
  tar -xvf gap2seq-${GAP2SEQ_VER}-x64-linux.tar.bz2 && \
  mv gap2seq-${GAP2SEQ_VER}-x64-linux/* /usr/local/bin

RUN Gap2Seq -h

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 \
    && tar xjvf samtools-${SAMTOOLS_VER}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VER} \
    && ./configure \
    && make \
    && make install

CMD ["bash"]
