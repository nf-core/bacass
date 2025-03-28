process NANOSTATS {

    tag "ASSAMBLE FLY LONG READS ${sample_code}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/675720c3645bd69ae8bdd55a601e38caa4bbd08f4db1567ab52da778b185c351/data' :
        'community.wave.seqera.io/library/pip_nanostat:69dc9b3e7f904176' }"


    input:

    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*"), emit: info_cov


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    # Generate stats QC
    NanoStat --fastq ${reads} > dedup_${prefix}_NanoStat.log

    """
}
