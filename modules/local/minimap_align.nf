// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/minimap2:2.21--h5bf99c6_0"
    }

    input:
    tuple val(meta), val(reads), file(longreads), file('reference')

    output:
    tuple val(meta), val(reads), file(longreads), path("*.paf"), emit: paf
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    minimap2 \\
        $options.args \\
        -t $task.cpus \\
        reference \\
        $longreads \\
        > ${prefix}.paf

    echo \$(minimap2 --version 2>&1) > ${software}.version.txt
    """
}
