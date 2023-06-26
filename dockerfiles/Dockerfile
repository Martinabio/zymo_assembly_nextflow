FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .

RUN /opt/conda/bin/mamba env create -f /environment.yml

ENV PATH=/opt/conda/envs/zymo_assembly_nextflow/bin:$PATH

CMD ["/bin/bash"]
