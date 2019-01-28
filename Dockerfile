FROM nfcore/base
LABEL authors="Andreas Wilm" \
      description="Docker image containing all requirements for nf-core/bacass pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-bacass-latest/bin:$PATH
