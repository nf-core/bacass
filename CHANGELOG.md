# nf-core/bacass: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0dev nf-core/bacass: "Aqua Platinum Zebrafish"

### `Changed`

### `Added`

- [#104](https://github.com/nf-core/bacass/pull/104), [#106](https://github.com/nf-core/bacass/pull/106) - Added dragonflye module and enbled draft genome polishing with short-reads.

### `Fixed`

### `Dependencies`

### `Deprecated`

## v2.1.0 nf-core/bacass: "Navy Steel Swordfish" - 2023/10/20

This version merges the nf-core template updates of v2.9 and v2.10, and updates modules or dependencies to ensure compatibility with the new template. Additionally, new modules have been added to process short-reads and perform gene annotation with Bakta.

### `Changed`

- [#86](https://github.com/nf-core/bacass/pull/86) - Update nf-core/bacass to the new nf-core 2.9 `TEMPLATE`.
- [#61](https://github.com/nf-core/bacass/issues/61) - Update local/modules to nf-core/modules (detailed below).
- [#91](https://github.com/nf-core/bacass/pull/91) - Update nf-core/bacass to the new nf-core 2.10 `TEMPLATE`.
- [#95](https://github.com/nf-core/bacass/pull/95) - Update MultiQC module to v1.17.

### `Added`

- [#86](https://github.com/nf-core/bacass/pull/86) - Added nf-core subworkflow for trimming and QC of short-reads [nf-core/fastq_trim_fastp_fastqc](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/fastq_trim_fastp_fastqc).
- [#88](https://github.com/nf-core/bacass/pull/88) - Added nf-validation on samplesheet
- [#93](https://github.com/nf-core/bacass/pull/93) - Added missing modules output to MultiQC. ( Fastp, PycoQC, Porechop, Quast, Kraken2, and Prokka).
- [#95](https://github.com/nf-core/bacass/pull/95) - Added subworkflow for gene annotation with Bakta.

### `Fixed`

- Fixed modules
  - Medaka: Medaka last version (see version update below) doesn't allow gzip compressed files. Add bgzip compression instead.
  - Dfast: fix overwriting issues detected when copying sample files from `work/` to `results/`

### `Dependencies`

- [#61](https://github.com/nf-core/bacass/issues/61) - Update local/modules to nf-core/modules plus version update.

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| Canu     | 2.1.1            | 2.2         |
| Minimap2 | 2.21             | 2.2         |
| Miniasm  | 0.3              | -           |
| Racon    | 1.4.20-1         | -           |

- Update already nf-core modules

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| Fastqc   | 0.11.9           | -           |
| Samtools | 1.13             | 2.1.2       |
| Kraken2  | 2.1.1            | 2.1.2       |
| Quast    | 5.0.2            | 5.2.0       |
| Prokka   | 1.14.6           | -           |
| Multiqc  | 1.10.1           | 1.15        |

- Refactor `local/modules` making them follow nf-core v2.9 structure/fashion.

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| Dfast      | 1.2.14           | -           |
| Medaka     | 1.4.3-0          | -           |
| Nanoplot   | 1.38.0           | 1.41.6      |
| Nanopolish | 0.13.2-5         | 0.14.0      |
| Pycoqc     | 2.5.2            | -           |
| Unicycler  | 0.4.8            | -           |

### `Deprecated`

- [#86](https://github.com/nf-core/bacass/pull/86) Replace depecated modules with nf-core/modules.

  - Replace `local/get_software_versions.nf` with `nf-core/custom/dumpsoftwareversions.nf`
  - Replace `local/skewer` by `nf-core/fastp` and wrap fastqc plus fastp into `subworkflows/nf-core/fastq_trim_fastp_fastqc`

## v2.0.0 nf-core/bacass: "Navy Steel Swordfish" 2021/08/27

### `Changed`

- [#56](https://github.com/nf-core/bacass/pull/56) - Switched to DSL2 & update to new nf-core 2.1 `TEMPLATE`
- [#56](https://github.com/nf-core/bacass/pull/56) - `--krakendb` now expects a `.tar.gz`/`.tgz` (compressed tar archive) directly from `https://benlangmead.github.io/aws-indexes/k2` instead of an uncompressed folder.

### `Added`

- [#56](https://github.com/nf-core/bacass/pull/56) - Added full size test dataset, two Zetaproteobacteria sequenced with Illumina MiSeq Reagent Kit V2, PE250, 3 to 4 million read pairs.

### `Fixed`

- [#51](https://github.com/nf-core/bacass/issues/51) - Fixed Unicycler

### `Dependencies`

- [#56](https://github.com/nf-core/bacass/pull/56) - Updated a bunch of dependencies (unchanged: FastQC, Miniasm, Prokka, Porechop, QUAST)
  - Unicycler from 0.4.4 to 0.4.8
  - Kraken2 from 2.0.9beta to 2.1.1
  - MultiQC from 1.9 to 1.10.1
  - PYCOQC from 2.5.0.23 to 2.5.2
  - Samtools from 1.11 to 1.13
  - Canu from 2.0 to 2.1.1-2
  - dfast from 1.2.10 to 1.2.14
  - Medaka from 1.1.2 to 1.4.3-0
  - Minimap 2 from 2.17 to 2.21
  - Nanoplot from 1.32.1 to 1.38.0
  - Nanopolish from 0.13.2 to 0.13.2-5
  - Racon from 1.4.13 to 1.4.20-1
  - Skewer from 0.2.2 to 0.2.2-3

### `Deprecated`

## v1.1.1 nf-core/bacass: "Green Aluminium Shark" 2020/11/05

This is basically a maintenance update that includes template updates, fixed environments and some minor bugfixes.

- Merged in nf-core/tools template v 1.10.2
- Updated dependencies
  - fastqc=0.11.8, 0.11.9
  - multiqc=1.8, 1.9
  - kraken2=2.0.8_beta, 2.0.9beta
  - prokka=1.14.5, 1.14.6
  - nanopolish=0.11.2, 0.13.2
  - parallel=20191122, 20200922
  - racon=1.4.10, 1.4.13
  - canu=1.9, 2.0
  - samtools=1.9, 1.11
  - nanoplot=1.28.1, 1.32.1
  - pycoqc=2.5.0.3, 2.5.0.23
- Switched out containers for many tools to make DSLv2 transition easier (escape from dependency hell)

## v1.1.0 nf-core/bacass: "Green Aluminium Shark" 2019/12/13

- Added support for hybrid assembly using Nanopore and Illumina Short Reads
- Added methods for long-read Nanopore data
  - Nanopolish, for polishing of Nanopore data with Illumina reads
  - Medaka, as alternative assembly polishing method
  - PoreChop, for quality trimming of Nanopore data
  - Nanoplot, for plotting quality metrics of Nanopore data
  - PycoQC, to QC Nanopore data
- Added multiple tools to assemble long-reads
  - Miniasm + Racon
  - Canu Assembler
  - Unicycler in Long read Mode
- Add alternative assembly annotation using DFAST
- Add social preview image

### Dependency updates

- Bumped Nextflow Version to 19.10.0

## Added tools

- DFAST
- PycoQC
- Nanoplot
- PoreChop
- Nanopolish

## v1.0.0 nf-core/bacass: "Green Tin Ant"

Initial release of nf-core/bacass, created with the [nf-core](http://nf-co.re/) template.

This pipeline is for bacterial assembly of next-generation sequencing reads. It can be used to quality trim your reads using [Skewer](https://github.com/relipmoc/skewer) and performs basic QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Afterwards, the pipeline performs read assembly using [Unicycler](https://github.com/rrwick/Unicycler) and assesses assembly quality using [QUAST](http://bioinf.spbau.ru/quast). Contamination of the assembly is checked using [Kraken2](https://ccb.jhu.edu/software/kraken2/) to verify sample purity. The resulting bacterial assembly is annotated using [Prokka](https://github.com/tseemann/prokka).

Furthermore, the pipeline creates various reports in the `results` directory specified, including a [MultiQC](https://multiqc.info) report summarizing some of the findings and software versions.
