# Containers

## `nfcore/bacass` container

Our main container is designed using [Conda](https://conda.io/) to install all tools used in nf-core/bacass.

## `Porechop` container

This container contains the porechop tool that is incompatible to be distributed with the standard container itself. Therefore, we have a secondary container that is used instead to run porechop only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.

## `QUAST` container

This container contains the QUAST tool that is incompatible to be distributed with the standard container itself. Therefore, we have a secondary container that is used instead to run QUAST only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.

## `DFAST` container

This container contains the DFAST tool for annotation of the assembled genome and cannot be distributed along with the standard container due to incompatible Python requirements with Prokka. Therefore, we have a secondary container that is used instead to run DFAST only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.

## `Medaka` container

This container contains the Medaka tool for polishing of the assembled genome and cannot be distributed along with the standard container due to incompatible Python requirements with PycoQC. Therefore, we have a secondary container that is used instead to run Medaka only. We utilize the Biocontainers container that is built using the bioconda recipe for this instead.
