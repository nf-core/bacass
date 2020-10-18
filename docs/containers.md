# Containers

## `nfcore/bacass` container

Some of the tools we use in bacass are inherently relying on software dependencies that are not installable in the same environment due to dependency conflicts. In these cases we utilize biocontainers. With a potential switch to Nextflow DSLv2 and therefore nextflow modules, we intend to have separated containers for all processes. This subsequently also means, that we cannot provide a single conda environment for the entire pipeline - please use Singularity or Docker to run your analysis - there is _NO_ conda support with the pipeline.
