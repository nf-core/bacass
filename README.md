# nf-core/bacass: Simple bacterial assembly and annotation pipeline

[![Build Status](https://travis-ci.com/nf-core/bacass.svg?branch=master)](https://travis-ci.com/nf-core/bacass)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A518.10.1-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bacass.svg)](https://hub.docker.com/r/nfcore/bacass)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

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

This pipeline is for bacterial assembly of next-gen sequencing reads. It quality-trims your reads, runs QC, including a
metagenomics classifier on your reads to check for purity, assembles the reads with
Unicycler, annotates the assembly with Prokka and runs Quast on the assembly for quality check.

## Credits

nf-core/bacass was originally written by Andreas Wilm.
