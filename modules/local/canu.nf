// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CANU {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    label 'process_high_memory'
    label 'error_retry'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'canu=2.1.1-2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/canu:2.1.1--h1b792b2_2"
    } else {
        container "quay.io/biocontainers/canu:2.1.1--h1b792b2_2"
    }

    input:
    tuple val(meta), val(reads), file(longreads)

    output:
    tuple val(meta), path('*_assembly.fasta') , emit: assembly
    tuple val(meta), path('*_assembly.report'), emit: log
    path  '*.version.txt'                     , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def genomeSize  = meta.genome_size == 'NA' ? "5m" : "${meta.genome_size}"
    """
    canu -p assembly -d canu_out \
        ${options.args} \
        genomeSize="${genomeSize}" -nanopore "${longreads}" \
        maxThreads="${task.cpus}" merylMemory="${task.memory.toGiga()}G" \
        merylThreads="${task.cpus}" hapThreads="${task.cpus}" batMemory="${task.memory.toGiga()}G" \
        redMemory="${task.memory.toGiga()}G" redThreads="${task.cpus}" \
        oeaMemory="${task.memory.toGiga()}G" oeaThreads="${task.cpus}" \
        corMemory="${task.memory.toGiga()}G" corThreads="${task.cpus}"
    mv canu_out/assembly.contigs.fasta ${prefix}_assembly.fasta
    mv canu_out/assembly.report ${prefix}_assembly.report

    echo \$(canu --version 2>&1) | sed -e 's/Canu //g' > ${software}.version.txt
    """
}
