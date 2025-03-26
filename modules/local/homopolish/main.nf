process HOMOPOLISH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homopolish:0.4.1--pyhdfd78af_1' :
        'biocontainers/homopolish:0.4.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(medaka_genome)
    tuple val(meta_gunzip), path(bacteria_sketch)
    
    output:
    tuple val(meta), path('*_genome_homopolished.fasta') , emit: assembly
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    homopolish polish \
        -a $medaka_genome \
        -s $bacteria_sketch \
        -m $params.homopolish_model \
        -o .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homopolish: \$( homopolish --version 2>&1 | sed 's/Homopolish VERSION: *//g' )
    END_VERSIONS

    """
}
