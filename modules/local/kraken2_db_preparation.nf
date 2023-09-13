process KRAKEN2_DB_PREPARATION {
    tag "${db.simpleName}"
    label 'process_low'

    conda 'conda-forge::sed=4.7'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path db

    output:
    tuple val("${db.simpleName}"), path("database"), emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/
    """
}
