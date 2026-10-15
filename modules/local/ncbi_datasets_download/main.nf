process NCBI_DATASETS_DOWNLOAD {
    tag "${accession}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/ncbi-datasets-cli_unzip:91afc97698fc876a' :
        'community.wave.seqera.io/library/ncbi-datasets-cli_unzip:91afc97698fc876a' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fna"), emit: fna
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    datasets download genome accession ${accession} \\
        --include genome,gff3 \\
        --filename dataset.zip \\
        $args

    ## Extract efficiently
    unzip -j dataset.zip "ncbi_dataset/data/${accession}*/*_genomic.fna" && mv *_genomic.fna ${accession}_genomic.fna
    unzip -j dataset.zip "ncbi_dataset/data/${accession}*/genomic.gff" && mv genomic.gff ${accession}.gff

    rm -f dataset.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: \$(datasets --version 2>&1 | grep -oP 'datasets version \\K[0-9.]+')
    END_VERSIONS
    """

    stub:
    """
    touch ${accession}_genomic.fna
    touch ${accession}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: \$(datasets --version 2>&1 | grep -oP 'datasets version \\K[0-9.]+')
    END_VERSIONS
    """
}
