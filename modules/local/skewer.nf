// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process SKEWER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "skewer=0.2.2-3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/skewer:0.2.2--hc9558a2_3"
    } else {
        container "quay.io/biocontainers/skewer:0.2.2--hc9558a2_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trm-cmb.R{1,2}.fastq.gz")   , emit: reads
    path("*.log"), emit: log
    //path "*.version.txt"                       , emit: version // TODO: output version!

    script:
    def software    = getSoftwareName(task.process)
    """
    # loop over readunits in pairs per sample
    pairno=0
    echo "${reads[0]} ${reads[1]}" | xargs -n2 | while read fq1 fq2; do
        skewer $options.args -t ${task.cpus} \$fq1 \$fq2;
    done

    # gzip, because skewer's -z returns an error
    gzip *.fastq

    cat \$(ls *trimmed-pair1.fastq.gz | sort) >> ${meta.id}_trm-cmb.R1.fastq.gz
    cat \$(ls *trimmed-pair2.fastq.gz | sort) >> ${meta.id}_trm-cmb.R2.fastq.gz

    #skewer --version
    """
}
