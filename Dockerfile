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
RUN conda --version

# install docker
RUN apt update -y
RUN apt install -y curl
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-latest.tgz | tar xvz -C /tmp/ && mv /tmp/docker/docker /usr/bin/docker


# install pecat
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN conda create -n pecat-env pecat=0.0.2 minimap2 racon perl samtools

#RUN conda init bash && . /root/.bashrc && conda activate pecat-env

ENV PATH /root/miniconda3/envs/pecat-env/bin:$PATH

WORKDIR /mnt
