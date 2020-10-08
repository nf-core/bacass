# ![nf-core/bacass](docs/images/nfcore-bacass_logo.png)

A simple bacterial assembly and annotation pipeline

[![GitHub Actions](https://github.com/nf-core/bacass/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/bacass/workflows/nf-core%20CI/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bacass.svg)](https://hub.docker.com/r/nfcore/bacass)
[![DOI](https://zenodo.org/badge/168486714.svg)](https://zenodo.org/badge/latestdoi/168486714)

## Introduction

### Short Read Assembly

This pipeline is primarily for bacterial assembly of next-generation sequencing reads. It can be used to quality trim your reads using [Skewer](https://github.com/relipmoc/skewer) and performs basic sequencing QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Afterwards, the pipeline performs read assembly using [Unicycler](https://github.com/rrwick/Unicycler). Contamination of the assembly is checked using [Kraken2](https://ccb.jhu.edu/software/kraken2/) to verify sample purity.

### Long Read Assembly

For users that only have Nanopore data, the pipeline quality trims these using [PoreChop](https://github.com/rrwick/Porechop) and assesses basic sequencing QC utilizing [NanoPlot](https://github.com/wdecoster/NanoPlot) and [PycoQC](https://github.com/a-slide/pycoQC).
The pipeline can then perform long read assembly utilizing [Unicycler](https://github.com/rrwick/Unicycler), [Miniasm](https://github.com/lh3/miniasm) in combination with [Racon](https://github.com/isovic/racon) or [Canu](https://github.com/marbl/canu). Long reads can be polished using specified Fast5 files with [NanoPolish](https://github.com/jts/nanopolish).

### Hybrid Assembly

For users specifying both short read and long read (NanoPore) data, the pipeline can perform a hybrid assembly approach utilizing [Unicycler](https://github.com/rrwick/Unicycler), taking the full set of information from short reads and long reads into account.

### Shared QC across all forms of assembly

In all cases, the assembly is assessed using [QUAST](http://bioinf.spbau.ru/quast).
The resulting bacterial assembly is furthermore annotated using [Prokka](https://github.com/tseemann/prokka).

In addition, the pipeline creates various reports in the `results` directory specified, including a [MultiQC](https://multiqc.info) report summarizing some of the findings and software versions.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple computing infrastructures in a portable manner. It comes with docker or singularity containers as well as conda environments, making installation trivial and results highly reproducible.

## Documentation

The nf-core/bacass pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

nf-core/bacass was originally written by [Andreas Wilm](https://github.com/andreas-wilm) and is currently maintained by [Alexander Peltzer](https://github.com/apeltzer).

## Citation

If you use nf-core/bacass for your analysis, please cite it using the following doi: [10.5281/zenodo.3574476](https://zenodo.org/badge/latestdoi/168486714)

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
