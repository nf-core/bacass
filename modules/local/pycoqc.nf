// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PYCOQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pycoqc=2.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pycoqc:2.5.2--py_0"
    } else {
        container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    }

    input:
    tuple val(meta), path(fast5)

    output:
    tuple val(meta), path('sequencing_summary.txt'), emit: summary
    path "*.html"        , emit: html
    path "*.json"        , emit: json
    path  "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    //Find out whether the sequencing_summary already exists
    if(file("${fast5}/sequencing_summary.txt").exists()){
        run_summary = "cp ${fast5}/sequencing_summary.txt ./sequencing_summary.txt"
    } else {
        run_summary =  "Fast5_to_seq_summary -f $fast5 -t ${task.cpus} -s './sequencing_summary.txt' --verbose_level 2"
    }
    //Barcodes available?
    barcode_me = file("${fast5}/barcoding_sequencing.txt").exists() ? "-b ${fast5}/barcoding_sequencing.txt" : ''
    """
    $run_summary

    pycoQC \\
        $options.args \\
        -f "sequencing_summary.txt" \\
        $barcode_me \\
        -o ${meta.id}_pycoqc.html \\
        -j ${meta.id}_pycoqc.json

    echo \$(pycoQC --version 2>&1) | sed 's/^.*pycoQC v//; s/ .*\$//' > ${software}.version.txt
    """
}
