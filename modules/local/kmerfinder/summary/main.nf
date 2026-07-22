process KMERFINDER_SUMMARY {
    tag "kmerfinder_summary"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path(report, stageAs: 'reports/*')

    output:
    path "*.csv"        , emit: summary
    path "*.yaml"       , emit: yaml
    path "versions.yml" , emit: versions

    script:
    """
    ## summarizing kmerfinder results
    kmerfinder_summary.py --path reports/ --output_bn kmerfinder.bn --output_csv kmerfinder_summary.csv

    ## Create a yaml file from csv
    csv_to_yaml.py -i kmerfinder_summary.csv -k 'sample_name' -op kmerfinder_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_CSV > kmerfinder_summary.csv
    sample_name,species,accession
    stub,Stub species,GCF_000000000.1_ASM000000v1
    END_CSV

    cat <<-END_YAML > kmerfinder_summary.yaml
    stub:
      species: Stub species
      accession: GCF_000000000.1_ASM000000v1
    END_YAML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
