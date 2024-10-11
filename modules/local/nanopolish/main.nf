process NANOPOLISH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--h773013f_3' :
        'biocontainers/nanopolish:0.14.0--h773013f_3' }"

    input:
    tuple val(meta), val(reads), file(longreads), file(assembly), file(bam), file(bai), file(fast5)

    output:
    tuple val(meta), file('polished_genome.fa') , emit: assembly
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    nanopolish index -d "${fast5}" "${longreads}"

    nanopolish variants \
        --consensus \
        -o polished.vcf \
        -r "${longreads}" \
        -b "${bam}" \
        -g "${assembly}" \
        -t "${task.cpus}" \
        --min-candidate-frequency 0.1 \
        $args

    nanopolish vcf2fasta -g "${assembly}" polished.vcf > polished_genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$( nanopolish --version | sed -e "s/nanopolish version //g" | head -n 1 )
    END_VERSIONS
    """
}
