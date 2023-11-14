process KMERFINDER_SUMMARY {
    tag "kmerfinder_summary"
    label 'process_low'

    conda "bioconda::python=3.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    path(reports, stageAs: 'reports/*')

    output:
    path "kmerfinder.csv"   , emit: summary
    path "versions.yml"     , emit: versions

    script:
    """
    kmerfinder_summary.py --path reports/ --output_bn kmerfinder.bn --output_csv kmerfinder.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
