// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NANOPOLISH {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    label 'process_high_memory'
    label 'error_retry'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'nanopolish=0.13.2-5' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanopolish:0.13.2--h8cec615_5"
    } else {
        container "quay.io/biocontainers/nanopolish:0.13.2--h8cec615_5"
    }

    input:
    tuple val(meta), val(reads), file(longreads), file(assembly), file(bai), file(fast5)

    output:
    tuple val(meta), file('polished_genome.fa'), emit: assembly
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    nanopolish index -d "${fast5}" "${longreads}"

    nanopolish variants \
        --consensus \
        -o polished.vcf \
        -r "${longreads}" \
        -b "${prefix}.bam" \
        -g "${assembly}" \
        -t "${task.cpus}" \
        --min-candidate-frequency 0.1

    nanopolish vcf2fasta -g "${assembly}" polished.vcf > polished_genome.fa

    nanopolish --version | sed -e "s/nanopolish version //g" | head -n 1 > ${software}.version.txt
    """
}
