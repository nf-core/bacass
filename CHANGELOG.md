# nf-core/bacass: Changelog

## v1.1.1 nf-core/bacass: "Green Aluminium Shark" 2020/11/05

This is basically a maintenance update that includes template updates, fixed environments and some minor bugfixes.

* Merged in nf-core/tools template v 1.10.2
* Updated dependencies
  * fastqc=0.11.8, 0.11.9
  * multiqc=1.8, 1.9
  * kraken2=2.0.8_beta, 2.0.9beta
  * prokka=1.14.5, 1.14.6
  * nanopolish=0.11.2, 0.13.2
  * parallel=20191122, 20200922
  * racon=1.4.10, 1.4.13
  * canu=1.9, 2.0
  * samtools=1.9, 1.11
  * nanoplot=1.28.1, 1.32.1
  * pycoqc=2.5.0.3, 2.5.0.23
* Switched out containers for many tools to make DSLv2 transition easier (escape from dependency hell)

## v1.1.0 nf-core/bacass: "Green Aluminium Shark" 2019/12/13

* Added support for hybrid assembly using Nanopore and Illumina Short Reads
* Added methods for long-read Nanopore data
  * Nanopolish, for polishing of Nanopore data with Illumina reads
  * Medaka, as alternative assembly polishing method
  * PoreChop, for quality trimming of Nanopore data
  * Nanoplot, for plotting quality metrics of Nanopore data
  * PycoQC, to QC Nanopore data
* Added multiple tools to assemble long-reads
  * Miniasm + Racon
  * Canu Assembler
  * Unicycler in Long read Mode
* Add alternative assembly annotation using DFAST
* Add social preview image

### Dependency updates

* Bumped Nextflow Version to 19.10.0

## Added tools

* DFAST
* PycoQC
* Nanoplot
* PoreChop
* Nanopolish

## v1.0.0 nf-core/bacass: "Green Tin Ant"

Initial release of nf-core/bacass, created with the [nf-core](http://nf-co.re/) template.

This pipeline is for bacterial assembly of next-generation sequencing reads. It can be used to quality trim your reads using [Skewer](https://github.com/relipmoc/skewer) and performs basic QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Afterwards, the pipeline performs read assembly using [Unicycler](https://github.com/rrwick/Unicycler) and assesses assembly quality using [QUAST](http://bioinf.spbau.ru/quast). Contamination of the assembly is checked using [Kraken2](https://ccb.jhu.edu/software/kraken2/) to verify sample purity. The resulting bacterial assembly is annotated using [Prokka](https://github.com/tseemann/prokka).

Furthermore, the pipeline creates various reports in the `results` directory specified, including a [MultiQC](https://multiqc.info) report summarizing some of the findings and software versions.
