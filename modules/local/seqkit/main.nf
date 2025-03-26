process SEQUIT {

    tag "$meta"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/902e5f138e9145c41c3e7848a9d8b7b7f8a21335ac76ac811b96cabcb3d277ad/data' :
        'community.wave.seqera.io/library/seqkit:2.10.0--03b4774218b4b7ef' }"

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("dedup_${meta.id}.fastq"), emit: seqkit_fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def valid_modes = ['-s', '-S', '-n', '-N']
    if (!valid_modes.contains(mode)) {
        throw new IllegalArgumentException("Unrecognised mode to run seqkit rmdup. Options: ${valid_modes.join(', ')}")
    }

    """
    seqkit rmdup ${mode} -i ${reads} -o dedup_${meta.id}.fastq
    """
}
