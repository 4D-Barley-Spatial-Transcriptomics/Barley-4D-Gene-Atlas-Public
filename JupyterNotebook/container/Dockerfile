# Based on the jupyter/r-notebook image
# https://github.com/jupyter/docker-stacks/blob/main/images/r-notebook/Dockerfile
ARG REGISTRY=quay.io
ARG OWNER=jupyter
ARG BASE_CONTAINER=$REGISTRY/$OWNER/r-notebook:x86_64-r-4.3.2
FROM $BASE_CONTAINER

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

# Seurat pre-requisites
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends libhdf5-dev libfontconfig1-dev \
    figlet \
    git-lfs \
    moreutils \
    vmtouch \
    locales \
    curl \
    sudo \
    nano \
    tree \
    build-essential \
    zip unzip wget \
    software-properties-common \
    libglpk40 \
    libfftw3-dev 

USER ${NB_UID}

WORKDIR "/tmp/work"

RUN mamba install --yes libcurl pkg-config r-igraph r-hdf5r r-fftw r-r.utils r-rfast

RUN mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "install.packages('Rfast2', repos = c(CRAN = 'https://cloud.r-project.org'))"

COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY=renv/library

RUN R -e "renv::restore()"

RUN R -e "options(repos = c(CRAN = 'https://cran.rstudio.com')); \
          install.packages('BiocManager'); \
          BiocManager::install('EBImage')"

