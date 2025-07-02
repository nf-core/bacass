process KMERFINDER_DOWNLOAD_REFERENCE {
    tag "${task.process}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0' :
        'biocontainers/requests:2.26.0' }"

    input:
    tuple val(refmeta), path(reports,  stageAs: 'reports/*')

    output:
    tuple val(refmeta), path("ncbi_dataset/data/*.fna") , emit: fna
    tuple val(refmeta), path("ncbi_dataset/data/*.gff") , emit: gff
    tuple val(refmeta), path("references_found.tsv")    , emit: references_tsv
    tuple val(refmeta), path("*.winner")                , emit: winner
    path "versions.yml"                                 , emit: versions

    script:
    """
    ## Find the common reference genome
    find_common_reference.py \\
        -d reports/ \\
        -o references_found.tsv

    ## Download the winner reference genome from the ncbi database
    download_reference.py \\
        -file references_found.tsv \\
        -out_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
