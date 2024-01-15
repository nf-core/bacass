process FIND_DOWNLOAD_REFERENCE {
    tag "${task.process}"
    label 'process_low'

    conda "conda-forge::requests=2.26.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0' :
        'biocontainers/requests:2.26.0' }"

    input:
    tuple val(meta), path(reports,  stageAs: 'reports/*')
    path(ncbi_metadata_db)

    output:
    tuple val(meta), path("*.fna.gz")              , emit: fna
    tuple val(meta), path("*.gff.gz")              , emit: gff
    tuple val(meta), path("*.faa.gz")              , emit: faa
    tuple val(meta), path("references_found.tsv")  , emit: references_tsv
    path "versions.yml"                            , emit: versions

    script:
    """
    ## Find the common reference genome
    find_common_reference.py \\
        -d reports/ \\
        -o references_found.tsv

    ## Download the winner reference genome from the ncbi database
    download_reference.py \\
        -file references_found.tsv \\
        -reference $ncbi_metadata_db \\
        -out_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
