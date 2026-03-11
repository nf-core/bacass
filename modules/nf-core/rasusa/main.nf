process RASUSA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rasusa:2.1.0--h4ac6f70_0' :
        'biocontainers/rasusa:2.1.0--h4ac6f70_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        if ("$reads" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
        """
        rasusa reads \\
            $args \\
            --output ${prefix}.fastq.gz \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rasusa: \$(rasusa --version | sed 's/rasusa //')
        END_VERSIONS
        """
    } else {
        if (reads[0].name == "${prefix}_1.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
        """
        rasusa reads \\
            $args \\
            --output ${prefix}_1.fastq.gz \\
            --output ${prefix}_2.fastq.gz \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rasusa: \$(rasusa --version | sed 's/rasusa //')
        END_VERSIONS
        """
    }

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def output_files = meta.single_end ? "${prefix}.fastq.gz" : "${prefix}_1.fastq.gz ${prefix}_2.fastq.gz"
    """
    touch $output_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version | sed 's/rasusa //')
    END_VERSIONS
    """
}
