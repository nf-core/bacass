process FLY {

    tag "ASSAMBLE FLY LONG READS ${meta}"
    
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa1c1e961de38d24cf36c424a8f4a9920ddd07b63fdb4cfa51c9e3a593c3c979/data' :
        'community.wave.seqera.io/library/flye:2.9.5--d577924c8416ccd8' }"
    
    input:

    tuple val(meta), path(reads)
    val mode
    val genomesize

    output:

    tuple val(meta),        path("*/assembly.fasta")        , emit: fly_assambly_tuple
    tuple val(meta),        path("*/assembly_graph.gfa")    , emit: grafic_assemble
    tuple val(meta),        path("*/assembly_info.txt")     , emit: info_cov
    path "versions.yml"                                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["-nanopore"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Fly. Options: ${valid_mode.join(', ')}" }

    """
    flye --nano-raw ${reads} --out-dir flye_output_${prefix} --genome-size ${genomesize} --plasmids --threads 16
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fly: \$(echo \$(fly --version 2>&1) | sed 's/^.*fly //; s/Using.*\$//' )
    END_VERSIONS
    """

}