# nf-core/bacass: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.5.0dev

### `Changed`

- [#221](https://github.com/nf-core/bacass/pull/221) update fastq_trim_fastp_fastqc subworkflow and its modules.

### `Added`

- [#222](https://github.com/nf-core/bacass/pull/222) Reinstalled Dragonflye 1.1.2.
- [#195](https://github.com/nf-core/bacass/pull/195) Update nf-core/bacass to the new nf-core 3.2.0 `TEMPLATE`.

### `Fixed`

- [#220](https://github.com/nf-core/bacass/pull/220) Fixed local environments.

### `Dependencies`

### `Deprecated`

## v2.4.0 nf-core/bacass: "Yellow Copper Crayfish" 2024/11/05

### `Changed`

- [#219](https://github.com/nf-core/bacass/pull/219) Required --assembly_type
- [#180](https://github.com/nf-core/bacass/pull/180) Bump version 2.4.0.
- [#169](https://github.com/nf-core/bacass/pull/169) Refactored long-reads polishing step.
- [#167](https://github.com/nf-core/bacass/pull/167) Remove params.save_merged as merged reads are not used in downstream analysis.
- [#159](https://github.com/nf-core/bacass/pull/159) Updated Kmerfinder module and increased memory.
- [#150](https://github.com/nf-core/bacass/pull/150) Replace local unicycler module with nf-core module + bump version.

### `Added`

- [#176](https://github.com/nf-core/bacass/pull/176) Update nf-core/bacass to nf-core-tools v3.0.2 `TEMPLATE`.
- [#166](https://github.com/nf-core/bacass/pull/166) Added FastQC after-trimming section to MultiQC report.
- [#158](https://github.com/nf-core/bacass/pull/158) Support automatic concatenation of FastQ files for the same sample.

### `Fixed`

- [#183](https://github.com/nf-core/bacass/pull/183) Fix DFAST issue in conda environment by updating its version.
- [#182](https://github.com/nf-core/bacass/pull/182) Uncommented required line to pass linting test in `--release` mode.
- [#179](https://github.com/nf-core/bacass/pull/179) Fixed matrix.test_name in linting and missing features from template 3.0.2.
- [#178](https://github.com/nf-core/bacass/pull/178) Fixed bakta running only for one sample.
- [#169](https://github.com/nf-core/bacass/pull/169) Fixed long reads polishing input channel.
- [#168](https://github.com/nf-core/bacass/pull/168) Fix wrong metadata in canu input channel.
- [#163](https://github.com/nf-core/bacass/pull/163) Fixed `params.save_merged` to properly save merged files.
- [#160](https://github.com/nf-core/bacass/pull/160) Fixed memory issues in KmerFinder, fixed handling of no species detected, and fixed handling of empty fasta files in the prokka/bakkta channel.
- [#157](https://github.com/nf-core/bacass/pull/157) Fixed corrupted zenodo URL of Kmerfinder database.
- [#154](https://github.com/nf-core/bacass/pull/154) Fixed kmerfinder script and increase resources to prevent memory issues.
- [#153](https://github.com/nf-core/bacass/pull/153) Update `.nf-core.yml` to fix files_unchanged section for accurate linting checks.

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| Dfast     | 1.2.20           | 1.3.2       |
| Unicycler | 0.4.8            | 0.5.0       |

### `Deprecated`

## v2.3.1 nf-core/bacass: "Navy Iron Oyster" 2024/06/24

### `Changed`

### `Added`

### `Fixed`

- [#147](https://github.com/nf-core/bacass/pull/147) Fixed input file errors related to samplesheets containing relative paths to symbolic links, addressing the 'not a valid path' error.

### `Dependencies`

### `Deprecated`

## v2.3.0 nf-core/bacass: "Navy Iron Oyster" 2024/06/12

### `Changed`

- [#135](https://github.com/nf-core/bacass/pull/135) Replaced nf-core MultiQC module with a custom MultiQC module.

### `Added`

- [#135](https://github.com/nf-core/bacass/pull/135) Implementation of KmerFinder subworkflow Custom Quast, and Custom MultiQC Reports:

  - Added KmerFinder subworkflow for read quality control, purity assessment, and sample grouping based on reference genome estimation.
  - Enhanced Quast Assembly QC to run both general and reference genome-based analyses when KmerFinder is invoked.
  - Implemented custom MultiQC module with multiqc_config.yml files for different assembly modes (short, long, hybrid).
  - Generated custom MultiQC HTML report consolidating metrics from KmerFinder, Quast, and other relevant sources.

- [#133](https://github.com/nf-core/bacass/pull/133) Update nf-core/bacass to the new nf-core 2.14.1 `TEMPLATE`.

### `Fixed`

- [#134](https://github.com/nf-core/bacass/pull/134) - Fixed samples reported of prokka/bakta in multiqc report.

- [#125](https://github.com/nf-core/bacass/pull/125) - Fixed conflicting settings in save_trimmed_fail parameter.

### `Dependencies`

### `Deprecated`

## v2.2.0 nf-core/bacass: "Aqua Platinum Zebrafish" 2024/03/30

### `Changed`

- [#111](https://github.com/nf-core/bacass/pull/111) - Update nf-core/bacass to 2.12, and [#117](https://github.com/nf-core/bacass/pull/117) to 2.13.1 `TEMPLATE`.

### `Added`

- [#104](https://github.com/nf-core/bacass/pull/104), [#106](https://github.com/nf-core/bacass/pull/106) - Added dragonflye module and enbled draft genome polishing with short-reads.

### `Fixed`

- [#123](https://github.com/nf-core/bacass/pull/123) - Add fallback to `download_pipeline.yml` when the pipeline does not support stub runs ([#2846](https://github.com/nf-core/tools/pull/2846))

### `Dependencies`

- [#120](https://github.com/nf-core/bacass/pull/120) - Update local and nf-core modules (version update an minnor code changes).

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| Bakta      | 1.8.2            | 1.9.3       |
| Canu       | 2.2              | -           |
| Dragonflye | 1.1.2            | -           |
| Fastp      | 0.23.4           | -           |
| Kraken2    | 2.1.2            | -           |
| Miniasm    | 0.3_r179         | -           |
| Minimap2   | 2.2              | 2.24        |
| Nanoplot   | 1.41.6           | -           |
| Porechop   | 0.2.4            | -           |
| Prokka     | 1.14.6           | -           |
| Quast      | 5.2.0            | -           |
| Racon      | 1.4.20           | -           |
| Samtools   | 1.17             | 1.19.2      |
| Untar      | 1.34             | -           |
| Dfast      | 1.2.20           | -           |
| Medaka     | 1.4.3-0          | -           |
| Nanopolish | 0.14.0           | -           |
| PycoQC     | 2.5.2            | -           |
| Unicycler  | 0.4.8            | -           |

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
