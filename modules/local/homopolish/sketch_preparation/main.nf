process HOMOPOLISH_SKETCH_PREPARATION {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/curl:7.80.0' :
        'biocontainers/curl:7.80.0' }"

    input:
    val(meta)
    path(url)

    output:
    tuple val(meta), path("bacteria.msh.gz"), emit: sketch
    path "versions.yml"                     , emit: versions

    script:
    """
    curl $params.homopolish_bacteria_sketch_url
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Homopolish_Sketch Bacteria: $params.homopolish_bacteria_last
    END_VERSIONS
    """
}
