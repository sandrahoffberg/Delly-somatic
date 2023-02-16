# hash:sha256:e5697a07d1358b759db59b9c8be71b602bb2558908fa8ea596c4186288590677
FROM registry.apps.codeocean.com/codeocean/jupyterlab:3.4.4-miniconda4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf=2.69-11.1 \
        automake=1:1.16.1-4ubuntu6 \
        awscli=1.18.69-1ubuntu0.20.04.1 \
        build-essential=12.8ubuntu1.1 \
        cmake=3.16.3-1ubuntu1.20.04.1 \
        git=1:2.25.1-1ubuntu3.10 \
        hdf5-tools=1.10.4+repack-11ubuntu1 \
        libboost-date-time-dev=1.71.0.0ubuntu2 \
        libboost-filesystem-dev=1.71.0.0ubuntu2 \
        libboost-iostreams-dev=1.71.0.0ubuntu2 \
        libboost-program-options-dev=1.71.0.0ubuntu2 \
        libboost-system-dev=1.71.0.0ubuntu2 \
        libbz2-dev=1.0.8-2 \
        libcurl4-gnutls-dev=7.68.0-1ubuntu2.15 \
        libgit2-dev=0.28.4+dfsg.1-2 \
        libhdf5-dev=1.10.4+repack-11ubuntu1 \
        liblzma-dev=5.2.4-1ubuntu1.1 \
        libncurses5-dev=6.2-0ubuntu2 \
        libssl-dev=1.1.1f-1ubuntu2.17 \
        pkg-config=0.29.1-0ubuntu4 \
        r-base=3.6.3-2 \
        wget=1.20.3-1ubuntu2 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.5 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        conda==4.12.0 \
    && conda clean -ya
RUN conda install -y --channel=conda-forge --channel=bioconda \
        bcftools==1.13 \
        htslib==1.13 \
        r-essentials==4.1 \
        r-ggplot2==3.4.0 \
        r-gtable==0.3.1 \
        r-reshape2==1.4.4 \
        r-scales==1.2.1 \
    && conda clean -ya

COPY postInstall /
RUN /postInstall
