process KMERFINDER_FIND_WINNER_REFERENCE {
    tag "${refmeta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(refmeta), path('reports/results_*.txt')

    output:
    tuple val(refmeta), path("references_found.tsv"), emit: references_tsv
    tuple val(refmeta), path("*.winner"), emit: winner
    tuple val(refmeta), val("${refmeta}_base_accession"), emit: base_accession  // ✅ Nuevo output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Find the common reference genome
    find_common_reference.py \\
        -d reports/ \\
        -o references_found.tsv

    ## Extract the winner accession from the TSV file
    FULL_ACCESSION=\$(head -n 1 references_found.tsv | grep -v '^#' | cut -f1)
    echo "Found full accession: \$FULL_ACCESSION"

    ## Extract base accession (remove assembly version part)
    ## GCF_002795805.1_ASM279580v1 → GCF_002795805.1
    BASE_ACCESSION=\$(echo "\$FULL_ACCESSION" | cut -d'_' -f1,2)
    echo "Base accession for datasets: \$BASE_ACCESSION"

    ## Create winner file with full accession (for compatibility)
    echo "\$FULL_ACCESSION" > \${FULL_ACCESSION}.winner

    ## Create base accession file for datasets tool
    echo "\$BASE_ACCESSION" > base_accession.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    ## Create stub files for testing
    echo "GCF_000000000.1_ASM000000v1\t1\tStub genome for testing" > references_found.tsv
    echo "GCF_000000000.1_ASM000000v1" > GCF_000000000.1_ASM000000v1.winner
    echo "GCF_000000000.1" > base_accession.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
    END_VERSIONS
    """
}
