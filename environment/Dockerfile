# hash:sha256:51799e397be28046f9e99c49605f7109c35ce914607305d325a4dae6e0263597
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf=2.71-2 \
        automake=1:1.16.5-1.3 \
        awscli=1.22.34-1 \
        build-essential=12.9ubuntu3 \
        cmake=3.22.1-1ubuntu1.22.04.1 \
        dos2unix=7.4.2-2 \
        git=1:2.34.1-1ubuntu1.9 \
        hdf5-tools=1.10.7+repack-4ubuntu2 \
        libboost-date-time-dev=1.74.0.3ubuntu7 \
        libboost-filesystem-dev=1.74.0.3ubuntu7 \
        libboost-iostreams-dev=1.74.0.3ubuntu7 \
        libboost-program-options-dev=1.74.0.3ubuntu7 \
        libboost-system-dev=1.74.0.3ubuntu7 \
        libbz2-dev=1.0.8-5build1 \
        libcurl4-gnutls-dev=7.81.0-1ubuntu1.10 \
        libgit2-dev=1.1.0+dfsg.1-4.1build1 \
        libhdf5-dev=1.10.7+repack-4ubuntu2 \
        liblzma-dev=5.2.5-2ubuntu1 \
        libncurses5-dev=6.3-2 \
        libssl-dev=3.0.2-0ubuntu1.9 \
        pkg-config=0.29.2-1ubuntu3 \
        r-base=4.1.2-1ubuntu2 \
        wget=1.21.2-2ubuntu1 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu9.2 \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && Rscript -e 'options(warn=2); install.packages("remotes")'

COPY postInstall /
RUN /postInstall
