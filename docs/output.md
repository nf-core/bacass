# nf-core/bacass: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Quality trimming and QC](#quality-trimming-and-qc)
    * [Short Read Trimming](#short-read-trimming)
    * [Short Read RAW QC](#short-read-raw-qc)
    * [Long Read Trimming](#long-read-trimming)
    * [Long Read RAW QC](#long-read-raw-qc)
* [Taxonomic classification](#taxonomic-classification)
* [Assembly Output](#assembly-output)
* [Assembly QC with QUAST](#assembly-qc-with-quast)
* [Annotation with Prokka](#annotation-with-prokka)
* [Annotation with DFAST](#annotation-with-dfast)
* [Report](#report)
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Quality trimming and QC

### Short Read Trimming

This step quality trims the end of reads, removes degenerate or too short reads and if needed,
combines reads coming from multiple sequencing runs.

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/trimming/shortreads/`
    * `*.fastq.gz`: Trimmed (and combined reads)

</details>

### Short Read RAW QC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/FastQC/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

![FastQC report](images/fastqc.png)

</details>

### Long Read Trimming

This step performs long read trimming on Nanopore input (if provided).

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/trimming/longreads/`
    * `trimmed.fastq.gz`: The trimmed FASTQ file

</details>

### Long Read RAW QC

These steps perform long read QC for input data (if provided).

Please refer to the documentation of [NanoPlot](https://github.com/wdecoster/NanoPlot) and [PycoQC](https://a-slide.github.io/pycoQC/) if you want to know more about the plots created by these tools.

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/QC_Longreads/`
    * `NanoPlot`: Various plots in HTML and PNG format
    * `PycoQC`

Example plot from Nanoplot:

![Nanoplot](images/nanoplot.png)

</details>

## Taxonomic classification

This QC step classifies your reads using [Kraken2](https://ccb.jhu.edu/software/kraken2/) a k-mer based approach. This helps to identify samples that have purity
issues. Ideally you will not want to assemble reads from samples that are contaminated or contain
multiple species. If you like to visualize the report, try
[Pavian](https://github.com/fbreitwieser/pavian) or [Krakey](http://krakey.info/).

<details markdown="1">
<summary>Output files</summary>

* `{sample}/Kraken2`
    * `{sample}.kraken2.report.txt`: Classification of short reads in the Kraken(1) report format.
    * `{sample}_longreads.kraken2.report.txt`: Classification of long reads in the Kraken(1) report format.

See [webpage](http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports) for more details.

Exemplary Kraken2 report screenshot:

![Kraken2 report](images/kraken2.png)

</details>

## Assembly Output

Trimmed reads are assembled with [Unicycler](https://github.com/rrwick/Unicycler) in `short` or `hybrid` assembly modes. For long-read assembly, there are also `canu` and `miniasm` available.
Unicycler is a pipeline on its own, which at least for Illumina reads mainly acts as a frontend to Spades with added polishing steps.

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/Unicycler`
    * `{sample}.scaffolds.fa`: Final assembly in fasta format
    * `{sample}.assembly.gfa`: Final assembly in Graphical Fragment Assembly (GFA) format
    * `{sample}.unicycler.log`: Log file summarizing steps and intermediate results on the Unicycler execution

Check out the [Unicycler documentation](https://github.com/rrwick/Unicycler) for more information on Unicycler output.

* `{sample_id}/Canu`
    * `{sample}_assembly.fasta`: Final assembly in fasta format
    * `{sample}_assembly.report`: Log file

Check out the [Canu documentation](https://canu.readthedocs.io/en/latest/index.html) for more information on Canu output.

* `{sample_id}/Miniasm`
    * `{sample}_assembly.fasta`: Assembly in fasta format
    * `{sample}_assembly_consensus.fasta`: Consensus assembly in fasta format (polished by Racon)

Check out the [Miniasm documentation](https://github.com/lh3/miniasm) for more information on Miniasm output.

</details>

## Assembly QC with QUAST

The assembly QC is performed with [QUAST](http://quast.sourceforge.net/quast) for all assemblies in one report. It reports multiple metrics including number of contigs, N50, lengths etc in form of an html report. It further creates an HTML file with integrated contig viewer (Icarus).

<details markdown="1">
<summary>Output files</summary>

* `QUAST`
    * `report.tsv`: QUAST's report in text format
* `QUAST/other_files`
    * `icarus.html`: QUAST's contig browser as HTML
    * `report.html`: QUAST assembly QC as HTML report
    * `report.pdf`: QUAST assembly QC as pdf

![QUAST QC](images/quast.png)

![Icarus](images/icarus.png)

</details>

## Annotation with Prokka

By default, the assembly is annotated with [Prokka](https://github.com/tseemann/prokka) which acts as frontend for several annotation tools and includes rRNA and ORF predictions.

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/Prokka/{sample_id}`
    * `{sample_id}.gff`: Annotation in gff format
    * `{sample_id}.txt`: Annotation in text format
    * `{sample_id}.faa`: Protein sequences in fasta format

See [Prokka's documentation](https://github.com/tseemann/prokka#output-files) for a full description of all output files.

![Prokka annotation](images/prokka.png)

</details>

## Annotation with DFAST

On request, the assembly is annotated with [DFAST](https://github.com/nigyta/dfast_core).

<details markdown="1">
<summary>Output files</summary>

* `{sample_id}/DFAST/RESULT_{dfast_profile_name}`
    * `genome.gff`: Annotation in gff format
    * `statistics.txt`: Annotation statistics in text format
    * `protein.faa`: Protein sequences in fasta format

</details>

## Report

Some pipeline results are visualised by [MultiQC](http://multiqc.info), which is a visualisation tool that generates a single HTML report summarising all samples in your project. Further statistics are available in within the report data directory.

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.tsv`.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>
