// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process PORECHOP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "porechop=0.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/porechop:0.2.4--py38hed8969a_1"
    } else {
        container "quay.io/biocontainers/porechop:0.2.4--py38hed8969a_1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('trimmed.fastq.gz'), emit: reads
    path "*.version.txt"                  , emit: version

    script:
    def software    = getSoftwareName(task.process)
    """
    porechop $options.args -i "${reads}" -t "${task.cpus}" -o trimmed.fastq.gz
    porechop --version > "${software}.version.txt"
    """
}
