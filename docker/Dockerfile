FROM mambaorg/micromamba:0.24.0 as app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

# metadata labels
LABEL base.image="mambaorg/micromamba:0.24.0"
LABEL dockerfile.version="1"
LABEL software="pbsc2"
LABEL software.version="1"
LABEL description="Conda environment for running PacBio's SARS-CoV-2 HiFi viral workflow"
LABEL website=""
LABEL license=""
LABEL license.url=""
LABEL maintainer1="Curtis Kapsak"
LABEL maintainer1.email="curtis.kapsak@theiagen.com"

# install dependencies; cleanup apt garbage; make /data directory
RUN apt-get update && apt-get install -y --no-install-recommends \
 wget \
 ca-certificates \
 git \
 procps && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/* && mkdir /data

# set the environment, can adjust PATH if needed. LC_ALL for singularity compatibility
ENV PATH="$PATH" \
 LC_ALL=C.UTF-8

# copy in conda environment yml file
COPY environment.yml /environment.yml

# install things into conda/mamba environment
RUN micromamba install -n base -f /environment.yml -y && \
 micromamba clean -a -y

# so that mamba/conda env is active when running below commands
ENV ENV_NAME="base"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

## INSTALL mimux - COPY IN BINARY?

# install VCFCons.py script manually
# RUN wget -O /usr/bin/VCFCons.py https://raw.githubusercontent.com/PacificBiosciences/CoSA/master/vcf/VCFCons.py && \
#  chmod +x /usr/bin/VCFCons.py

# install CoSA (which includes VCFCons.py)
RUN git clone https://github.com/Magdoll/CoSA.git && \
 cd CoSA && \
 python setup.py build && \
 python setup.py install

WORKDIR /data

# hardcode executables from base conda env into the PATH variable
ENV PATH="${PATH}:/opt/conda/bin/"

# new base for testing
FROM app as test

# so that mamba/conda env is active when running below commands
ENV ENV_NAME="base"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# simple tests to see versions and help options for important tools
RUN micromamba list && \
 samtools --help && \
 bcftools --help && \
 pbmm2 --help && \
 VCFCons.py --help
