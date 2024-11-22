FROM ubuntu:20.04

LABEL maintainer "epil02@uw.edu"

ARG MUMMER_VERSION=4.0.0rc1
ARG BIOPYTHON_VERSION=1.78
ARG MUMMER_VER="4.0.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=$PATH:/usr/local/bin

# Update package list and install necessary packages
RUN apt-get update && apt-get install -y --no-install-recommends \
  wget \
  git \
  libncurses5-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  zlib1g-dev \
  libssl-dev \
  gcc \
  make \
  perl \
  python3 \
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

# Install Biopython
RUN python3 -m pip install biopython==${BIOPYTHON_VERSION}

# Download and install MUMMER
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

CMD ["bash"]
