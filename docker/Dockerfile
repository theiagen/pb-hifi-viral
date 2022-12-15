FROM mambaorg/micromamba:1.1.0 as app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

# metadata labels
LABEL base.image="mambaorg/micromamba:1.1.0"
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

# copy in conda environment yml file
COPY environment.yml /environment.yml

# install things into base conda/mamba environment
RUN micromamba install -n base -f /environment.yml -y && \
 micromamba clean -a -y

# so that mamba/conda env is active when running below commands
ENV ENV_NAME="base"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# install VCFCons.py script manually
# RUN wget -O /usr/bin/VCFCons.py https://raw.githubusercontent.com/PacificBiosciences/CoSA/master/vcf/VCFCons.py && \
#  chmod +x /usr/bin/VCFCons.py

# install CoSA (which includes VCFCons.py)
RUN git clone https://github.com/Magdoll/CoSA.git && \
 cd CoSA && \
 python setup.py build && \
 python setup.py install

# copy in mimux binary
COPY mimux /usr/local/bin

# download WuHan-1 reference genome (MN908947.3)
RUN mkdir /reference-data && cd /reference-data/ && \
 wget -O MN908947.3.SARS-CoV-2.reference.fasta https://raw.githubusercontent.com/epi2me-labs/wf-artic/master/data/primer_schemes/SARS-CoV-2/ARTIC/V4.1/SARS-CoV-2.reference.fasta

WORKDIR /data

# hardcode executables from base conda env into the PATH variable; LC_ALL for singularity compatibility
ENV PATH="${PATH}:/opt/conda/bin/" \
 LC_ALL=C.UTF-8

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
 VCFCons.py --help && \
 mimux -h
