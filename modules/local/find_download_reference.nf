process FIND_DOWNLOAD_REFERENCE {
    tag "${task.process}"
    label 'process_low'

    conda "conda-forge::requests=2.26.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0' :
        'biocontainers/requests:2.26.0' }"

    input:
    tuple val(meta), path(reports,  stageAs: 'reports/*')
    path(ncbi_reference)

    output:
    path "references_found.tsv" , emit: target_references_tsv
    path "*.fna.gz"             , emit: fna
    path "*.gff.gz"             , emit: gff
    path "*.faa.gz"             , emit: faa
    path "versions.yml"         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.refseq}"
    """
    find_common_reference.py \\
        -d reports/ \\
        -o references_found.tsv

    download_reference.py \\
        -file references_found.tsv \\
        -reference $ncbi_reference \\
        -out_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
