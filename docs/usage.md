# nf-core/bacass: Usage

## Table of contents

* [nf-core/bacass: Usage](#nf-corebacass-usage)
  * [Table of contents](#table-of-contents)
  * [General Nextflow info](#general-nextflow-info)
  * [Running the pipeline](#running-the-pipeline)
    * [Updating the pipeline](#updating-the-pipeline)
    * [Reproducibility](#reproducibility)
  * [Main Nextflow arguments](#main-nextflow-arguments)
    * [`-profile`](#profile)
  * [Main Pipeline Arguments](#main-pipeline-arguments)
    * [`--annotation_tool`](#annotationtool)
    * [`--assembler`](#assembler)
    * [`--assembly_type`](#assemblytype)
    * [`--dfast_config`](#dfastconfig)
    * [`--input`](#input)
    * [`--kraken2db`](#kraken2db)
    * [`--prokka_args`](#prokkaargs)
    * [`--unicycler_args`](#unicyclerargs)
  * [Skipping Options](#skipping-options)
    * [`--skip_annotation`](#skipannotation)
    * [`--skip_kraken2`](#skipkraken2)
    * [`--skip_nanopolish`](#skipnanopolish)
    * [`--skip_pycoqc`](#skippycoqc)
  * [Job resources](#job-resources)
    * [Automatic resubmission](#automatic-resubmission)
    * [Custom resource requests](#custom-resource-requests)
  * [AWS Batch specific parameters](#aws-batch-specific-parameters)
    * [`--awsqueue`](#awsqueue)
    * [`--awsregion`](#awsregion)
  * [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#outdir)
    * [`--email`](#email)
    * [`-name`](#name)
    * [`-resume`](#resume)
    * [`-c`](#c)
    * [`--custom_config_version`](#customconfigversion)
    * [`--custom_config_base`](#customconfigbase)
    * [`--max_memory`](#maxmemory)
    * [`--max_time`](#maxtime)
    * [`--max_cpus`](#maxcpus)
    * [`--plaintext_email`](#plaintextemail)
    * [`--monochrome_logs`](#monochromelogs)
    * [`--multiqc_config`](#multiqcconfig)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/bacass --input design.tsv -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/bacass
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/bacass releases page](https://github.com/nf-core/bacass/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main Nextflow arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from DockerHub: [`nfcore/bacass`](http://hub.docker.com/r/nfcore/bacass/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

## Main Pipeline Arguments

### `--annotation_tool`

The annotation method to annotate the final assembly. Default choice is `prokka`, but the `dfast` tool is also available. For the latter, make sure to create your specific config if you're not happy with the default one provided. See [#dfast_config](#dfastconfig) to find out how.

### `--assembler`

The assembler to use for assembly. Available options are `Unicycler`, `Canu`, `Miniasm`. The latter two are only available for long-read data, whereas Unicycler can be used for short or hybrid assembly projects.

### `--assembly_type`

This adjusts the type of assembly done with the input data and can be any of `long`, `short` or `hybrid`. Short & Hybrid assembly will always run Unicycler, whereas long-read assembly can be configured separately using the `--assembler` parameter.

### `--dfast_config`

Specifies a configuration file for the [DFAST](https://github.com/nigyta/dfast_core) annotation method. This can be used instead of PROKKA if required to specify a specific config file for annotation. If you want to know how to create your config file, please refer to the [DFAST](https://github.com/nigyta/dfast_core) readme on how to create one.

### `--input`

Use this to specify the location of your input design file. For example:

```bash
--input 'design_hybrid.tsv'
```

An example of properly formatted input files can be found at the [nf-core/testData](https://github.com/nf-core/test-datasets/tree/bacass). Exemplarily, this is the input used for a hybrid assembly in testing:

```bash
ID R1 R2 LongFastQ Fast5 GenomeSize
ERR044595 https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_1.fastq.gz https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_2.fastq.gz https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz NA 2.8m
```

* `ID` The identifier to use for handling the dataset
* `R1` The forward reads in case of available short-read data
* `R2` The reverse reads in case of available short-read data
* `LongFastQ` The long*read FastQ file with reads in FASTQ format
* `Fast5` The folder containing the basecalled FAST5 files
* `GenomeSize` The expected genome size of the assembly. Only used by the canu assembler.

Missing values (e.g. FAST5 folder in case of short reads) can be omitted by using a `NA` in the TSV file. The pipeline will handle such cases appropriately then.

### `--kraken2db`

Path to Kraken2 database.
See [Kraken2 homepage](https://ccb.jhu.edu/software/kraken2/index.shtml#downloads) for download
links. Minikraken2 8GB is a reasonable choice, since we run Kraken here mainly just to check for
sample purity.

### `--prokka_args`

This advanced option allows you to pass extra arguments to Prokka (e.g. `" --rfam"` or `" --genus name"`). For this to work you need to quote the arguments and add at least one space

### `--unicycler_args`

This advanced option allows you to pass extra arguments to Unicycler (e.g. `"--mode conservative"` or `"--no_correct"`). For this to work you need to quote the arguments and add at least one space.

## Skipping Options

### `--skip_annotation`

Skip annotating the assembly with Prokka.

### `--skip_kraken2`

Skip running Kraken2 classifier on reads.

### `--skip_nanopolish`

Skip polishing the long-read assembly with FAST5 input. Will not affect short/hybrid assemblies.

### `--skip_pycoqc`

Skip running `PycoQC` on long read input.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [`Slack`](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
