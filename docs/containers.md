# Containers

## `nfcore/bacass` container

Our main container is designed using [Conda](https://conda.io/) to install all tools used in nf-core/bacass.

## `Porechop` container

This container contains the porechop tool that is incompatible to be distributed with the standard container itself. Therefore, we have a secondary container that is used instead to run porechop only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.

## `QUAST` container

This container contains the QUAST tool that is incompatible to be distributed with the standard container itself. Therefore, we have a secondary container that is used instead to run QUAST only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.
