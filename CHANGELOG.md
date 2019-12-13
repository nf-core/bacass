# nf-core/bacass: Changelog

## v1.1.0dev nf-core/bacass: "Green Aluminium Shark"

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
