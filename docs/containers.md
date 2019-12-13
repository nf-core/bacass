# Containers

## `nfcore/bacass` container

Our main container is designed using [Conda](https://conda.io/) to install all tools used in nf-core/bacass. All other dedicated containers contain individual tools that are not possible to be bundled with other software in the main container due to dependency collisions / incompatibilities. In this case, we use Biocontainers:

* Porechop
* QUAST
* DFAST
* Medaka
