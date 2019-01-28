From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Andreas Wilm
    DESCRIPTION Singularity image containing all requirements for the nf-core/bacass pipeline
    VERSION latest

%environment
    PATH=/opt/conda/envs/nf-core-bacass-latest/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
