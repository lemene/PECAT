FROM ubuntu:bionic

RUN apt-get update
RUN apt-get -y install curl
RUN apt-get install sudo

# Install Miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# install conda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# install docker
RUN apt update -y
RUN apt install -y curl
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-latest.tgz | tar xvz -C /tmp/ && mv /tmp/docker/docker /usr/bin/docker


# install pecat
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

ENV PATH /root/miniconda3/envs/pecat-env/bin:$PATH

RUN conda create -n pecat-env  minimap2=2.24 racon=1.5 perl=5.32 samtools=1.17 python=3.11

RUN wget https://github.com/lemene/PECAT/releases/download/v0.0.3/pecat_v0.0.3_d1e5be8.tar.gz
RUN mkdir -p /root/miniconda3/envs/pecat-env/share
RUN tar -zxf pecat_v0.0.3_d1e5be8.tar.gz -C /root/miniconda3/envs/pecat-env/share

ENV PATH /root/miniconda3/envs/pecat-env/share/pecat_v0.0.3/build/bin:$PATH

WORKDIR /mnt
