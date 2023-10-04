FROM continuumio/miniconda3:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# install docker
RUN apt update -y
RUN apt install -y curl
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-latest.tgz | tar xvz -C /tmp/ && mv /tmp/docker/docker /usr/bin/docker


# install third-party tools
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

RUN conda create -n pecat-env  minimap2=2.24 racon=1.5 perl=5.32 samtools=1.17 python=3.11 singularity=3.8
ENV PATH /opt/conda/envs/pecat-env/bin:$PATH

# install pecat
RUN export PECAT_BIN=pecat_v0.0.3_4661bd5.tar.gz && \
    wget https://github.com/lemene/PECAT/releases/download/v0.0.3/${PECAT_BIN} && \
    tar -zxf ${PECAT_BIN} -C /opt/conda/envs/pecat-env/share && \
    rm ${PECAT_BIN}
ENV PATH /opt/conda/envs/pecat-env/share/pecat_v0.0.3/build/bin:$PATH

WORKDIR /mnt
