// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIASM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::miniasm=0.3_r179' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/miniasm:0.3_r179--h5bf99c6_2"
    } else {
        container "quay.io/biocontainers/miniasm:0.3_r179--h5bf99c6_2"
    }

    input:
    tuple val(meta), val(reads), file(longreads), path(paf)

    output:
    tuple val(meta), path('*_assembly.fasta') , emit: assembly
    //tuple val(meta), path('*_assembly.report'), emit: log
    path  '*.version.txt'                     , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    miniasm -f "${longreads}" "${paf}" > "${longreads}.gfa"
    awk '/^S/{print ">"\$2"\\n"\$3}' "${longreads}.gfa" | fold > ${prefix}_assembly.fasta

    echo \$(miniasm -V 2>&1) > ${software}.version.txt
    """
}
