FROM condaforge/mambaforge:latest AS conda

ENV DEBIAN_FRONTEND noninteractive 

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    zlib1g-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/GaetanBenoitDev/metaMDBG.git

RUN /opt/conda/bin/mamba install -n base -c bioconda -c defaults gsl=2.7=he838d99_0 wfmash=0.8.2 samtools=1.6 minimap2=2.24 && conda env config vars set CPATH=${CONDA_PREFIX}/include:${CPATH}

RUN conda clean -afy

RUN mkdir metaMDBG/build && cd metaMDBG/build && cmake .. && make -j 3

ENV PATH="$PATH:/metaMDBG/build/bin"

CMD ["/bin/bash"]