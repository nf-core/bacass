FROM nfcore/base
LABEL authors="Andreas Wilm, Alexander Peltzer" \
      description="Docker image containing Porechop for the nf-core/bacass pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# for bandage :/ otherwise it complains about missing libGL.so.1
RUN apt-get install -y libgl1-mesa-glx && apt-get clean -y
ENV PATH /opt/conda/envs/bacassporechop/bin:$PATH