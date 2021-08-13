// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEDAKA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'medaka=1.4.3-0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/medaka:1.4.3--py38h130def0_0"
    } else {
        container "quay.io/biocontainers/medaka:1.4.3--py38h130def0_0"
    }

    input:
    tuple val(meta), file(assembly), val(reads), file(longreads)

    output:
    tuple val(meta), path('*_polished_genome.fa'), emit: assembly
    path  '*.version.txt'                        , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    medaka_consensus ${options.args} \
        -i ${longreads} \
        -d ${assembly} \
        -o "${prefix}_polished_genome.fa" \
        -t ${task.cpus}

    echo \$(medaka --version 2>&1) | sed -e 's/medaka //g' > ${software}.version.txt
    """
}
