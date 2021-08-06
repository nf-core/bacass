// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DFAST {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "dfast=1.2.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/dfast:1.2.14--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/dfast:1.2.14--h2e03b76_0"
    }

    input:
    tuple val(meta), path(fasta)
    file (config)

    output:
    tuple val(meta), path("RESULT*"), emit: reads
    path "*.version.txt"            , emit: version

    script:
    def software    = getSoftwareName(task.process)
    """
    dfast_file_downloader.py --protein dfast --dbroot .
    dfast --genome ${fasta} --config $config
    dfast --version | sed -e "s/DFAST ver. //g"  > "${software}.version.txt"
    """
}
