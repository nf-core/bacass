FROM nfcore/base
LABEL authors="Andreas Wilm" \
      description="Docker image containing all requirements for nf-core/bacass pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# for bandage :/ otherwise it complains about missing libGL.so.1
RUN apt-get install -y libgl1-mesa-glx && apt-get clean -y
ENV PATH /opt/conda/envs/nf-core-bacass-latest/bin:$PATH
