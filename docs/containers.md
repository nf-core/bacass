# Containers

## `nfcore/bacass` container

Our main container is designed using [Conda](https://conda.io/) to install all tools used in nf-core/bacass. All other dedicated containers contain individual tools that are not possible to be bundled with other software in the main container due to dependency collisions / incompatibilities. In these cases we utilize biocontainers for some steps. With a potential switch to Nextflow DSLv2 and therefore nextflow modules, we intend to have separated containers per process. This subsequently also means, that we cannot provide a single conda environment for the entire pipeline - please use Singularity or Docker to run your analysis!
