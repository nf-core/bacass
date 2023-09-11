process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda 'medaka=1.4.3-0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.3--py38h130def0_0' :
        'biocontainers/medaka:1.4.3--py38h130def0_0' }"

    input:
    tuple val(meta), file(longreads), file(assembly)

    output:
    tuple val(meta), path('*_polished_genome.fa')   , emit: assembly
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def prefix                  = task.ext.prefix ?: "${meta.id}"
    def reads_bgzip_command     = ("$longreads".endsWith('.gz')) ? "zcat $longreads | bgzip -c > ${prefix}.fastq.bgz" : ''
    def assembly_bgzip_command  = ("$assembly".endsWith('.gz'))  ? "zcat $assembly  | bgzip -c > ${prefix}.fasta.bgz" : ''
    if ("$longreads".endsWith('.gz')) { reads_bgzip_out     = "${prefix}.fastq.bgz"} else { reads_bgzip_out    = null }
    if ("$assembly".endsWith('.gz'))  { assembly_bgzip_out  = "${prefix}.fasta.bgz"} else { assembly_bgzip_out = null }

    """
    # Recompress with bgzip
    $reads_bgzip_command
    $assembly_bgzip_command

    medaka_consensus $args \
        -i ${ reads_bgzip_out ?: longreads } \
        -d ${ assembly_bgzip_out ?: assembly } \
        -o "${prefix}_polished_genome.fa" \
        -t $task.cpus
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
