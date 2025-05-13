<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-bacass_logo_dark.png">
    <img alt="nf-core/bacass" src="docs/images/nf-core-bacass_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/bacass/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/bacass/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/bacass/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/bacass/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/bacass/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2669428-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2669428)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/bacass)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23bacass-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/bacass)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/bacass** is a bioinformatics best-practice analysis pipeline for simple bacterial assembly and annotation. The pipeline is able to assemble short reads, long reads, or a mixture of short and long reads (hybrid assembly).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/bacass/results).

## Pipeline summary

### Short Read Assembly

This pipeline is primarily for bacterial assembly of next-generation sequencing reads. It can be used to quality trim your reads using [FastP](https://github.com/OpenGene/fastp) and performs basic sequencing QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Afterwards, the pipeline performs read assembly using [Unicycler](https://github.com/rrwick/Unicycler). Contamination of the assembly is checked using [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Kmerfinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/) to verify sample purity.

### Long Read Assembly

For users that only have Nanopore data, the pipeline quality trims these using [PoreChop](https://github.com/rrwick/Porechop) or filter long reads by quality using [Filtlong](https://github.com/rrwick/Filtlong) and assesses basic sequencing QC utilizing [NanoPlot](https://github.com/wdecoster/NanoPlot) and [PycoQC](https://github.com/a-slide/pycoQC). Contamination of the assembly is checked using [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Kmerfinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/) to verify sample purity.

The pipeline can then perform long read assembly utilizing [Unicycler](https://github.com/rrwick/Unicycler), [Miniasm](https://github.com/lh3/miniasm) in combination with [Racon](https://github.com/isovic/racon), [Canu](https://github.com/marbl/canu) or [Flye](https://github.com/fenderglass/Flye) by using the [Dragonflye](https://github.com/rpetit3/dragonflye)(\*) pipeline. Long reads assembly can be polished using [Medaka](https://github.com/nanoporetech/medaka) or [NanoPolish](https://github.com/jts/nanopolish) with Fast5 files.

> [!NOTE]
> Dragonflye is a comprehensive pipeline designed for genome assembly of Oxford Nanopore Reads. It facilitates the utilization of Flye (default), Miniasm, and Raven assemblers, along with Racon (default) and Medaka polishers. For more information, visit the [Dragonflye GitHub](https://github.com/rpetit3/dragonflye) repository.

### Hybrid Assembly

For users specifying both short read and long read (NanoPore) data, the pipeline can perform a hybrid assembly approach utilizing [Unicycler](https://github.com/rrwick/Unicycler) (short read assembly followed by gap closing with long reads) or [Dragonflye](https://github.com/rpetit3/dragonflye) (long read assembly followed by polishing with short reads), taking the full set of information from short reads and long reads into account.

### Assembly QC and annotation

In all cases, the assembly is assessed using [QUAST](http://bioinf.spbau.ru/quast) and [BUSCO](https://busco.ezlab.org/). The resulting bacterial assembly is furthermore annotated using [Prokka](https://github.com/tseemann/prokka), [Bakta](https://github.com/oschwengers/bakta) or [DFAST](https://github.com/nigyta/dfast_core).

If Kmerfinder is invoked, the pipeline will group samples according to the [Kmerfinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/)-estimated reference genomes. Afterwards, two QUAST steps will be carried out: an initial ('general') [QUAST](http://bioinf.spbau.ru/quast) of all samples without reference genomes, and subsequently, a 'by reference genome' [QUAST](http://bioinf.spbau.ru/quast) to aggregate samples with their reference genomes.

> [!NOTE]
> This scenario is supported when [Kmerfinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/) analysis is performed only.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.tsv`:

```tsv
ID      R1                            R2                            LongFastQ                    Fast5    GenomeSize
shortreads      ./data/S1_R1.fastq.gz       ./data/S1_R2.fastq.gz       NA                            NA      NA
longreads       NA                          NA                          ./data/S1_long_fastq.gz      ./data/FAST5  2.8m
shortNlong      ./data/S1_R1.fastq.gz       ./data/S1_R2.fastq.gz       ./data/S1_long_fastq.gz      ./data/FAST5  2.8m

```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Short read assembly with Unicycler, `--kraken2db` can be any [compressed database (`.tar.gz`/`.tgz`)](https://benlangmead.github.io/aws-indexes/k2):

```console
nextflow run nf-core/bacass -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.tsv --assembly_type 'short' --kraken2db "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz"
```

Long read assembly with Miniasm:

```console
nextflow run nf-core/bacass -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.tsv --assembly_type 'long' --assembler 'miniasm' --kraken2db "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz"
```

```bash
nextflow run nf-core/bacass \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.tsv \
  --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/bacass/usage) and the [parameter documentation](https://nf-co.re/bacass/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/bacass/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/bacass/output).

## Credits

nf-core/bacass was initiated by [Andreas Wilm](https://github.com/andreas-wilm), originally written by [Alex Peltzer](https://github.com/apeltzer) (DSL1), rewritten by [Daniel Straub](https://github.com/d4straub) (DSL2) and maintained by [Daniel Valle-Millares](https://github.com/Daniel-VM).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#bacass` channel](https://nfcore.slack.com/channels/bacass) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/bacass for your analysis, please cite it using the following doi: [10.5281/zenodo.2669428](https://doi.org/10.5281/zenodo.2669428)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
