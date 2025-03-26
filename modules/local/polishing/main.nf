process POLISHING_PROCESS {
    tag { "Polishing ${meta.id} - ${max_rounds} rounds" }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba89ca084feaa591520bf2f40c7f16238abfbb1e6bd70d4813fdb103a3514705/data':
        'community.wave.seqera.io/library/minimap2_racon:5f257adb6aaf9096' }"

    input:
    tuple val(meta), path(reads), path(assembly_fasta), val(max_rounds)

    output:
    tuple val(meta), path("final_polished_${meta.id}.fasta")

    script:
    """
    input_fasta=${assembly_fasta}

    for round in \$(seq 1 ${max_rounds}); do
        aln_file="aln_round\${round}_${meta.id}.sam"
        output_fasta="racon_round\${round}_${meta.id}.fasta"

        minimap2 -ax map-ont \${input_fasta} ${reads} -t 8 > \${aln_file}
        racon ${reads} \${aln_file} \${input_fasta} --threads 8 > \${output_fasta}

        if [ ! -f \${output_fasta} ]; then
            echo "Error: \${output_fasta} was not generated. Exiting." >&2
            exit 1
        fi

        input_fasta=\${output_fasta}
    done

    mv \${input_fasta} final_polished_${meta.id}.fasta
    """
}