# ![nf-core/bacass](docs/images/nfcore-bacass_logo.png)

A simple bacterial assembly and annotation pipeline

[![GitHub Actions CI Status](https://github.com/nf-core/bacass/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/bacass/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/bacass/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/bacass/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bacass.svg)](https://hub.docker.com/r/nfcore/bacass)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23bacass-4A154B?logo=slack)](https://nfcore.slack.com/channels/bacass)
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

The nf-core/bacass pipeline comes with documentation about the pipeline which you can read at [https://nf-core/bacass/docs](https://nf-core/bacass/docs) or find in the [`docs/` directory](docs).

## Credits

nf-core/bacass was originally written by Andreas Wilm, Alexander Peltzer.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#bacass` channel](https://nfcore.slack.com/channels/bacass) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use nf-core/bacass for your analysis, please cite it using the following doi: [10.5281/zenodo.3574476](https://zenodo.org/badge/latestdoi/168486714)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
