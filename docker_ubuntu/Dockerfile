FROM ubuntu:20.04

ARG PICARD_VERSION=3.30

ENV PICARD=/usr/local/bin/picard.jar
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \    

openjdk-11-jdk bc python3 wget libncurses5-dev bzip2 gcc g++ make zlib1g-dev libbz2-dev liblzma-dev unzip autoconf\

&& cd /tmp \
    && wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
    && mv picard.jar /usr/local/bin \

CMD ["bash"]
