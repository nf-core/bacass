process DFAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dfast:1.3.2--h43eeafb_0' :
        'biocontainers/dfast:1.3.2--h43eeafb_0' }"

    input:
    tuple val(meta), path(fasta)
    file (config)

    output:
    tuple val(meta), path("*_results"), emit: annotation
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    dfast_file_downloader.py \\
        $args \\
        --protein dfast \\
        --dbroot .

    dfast \\
        $args2 \\
        --genome ${fasta} \\
        --config $config

    mv RESULT_TEST/ ${prefix}_results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dfast: \$( dfast --version | sed -e "s/DFAST ver. //g" )
    END_VERSIONS
    """
}
