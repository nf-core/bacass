process CUSTOM_MULTIQC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path multiqc_logo
    path workflow_summary
    path methods_description
    path software_versions
    path ('fastqc/*')
    path ('fastqc_trim/*')
    path ('fastp/*')
    path ('nanoplot??/*')
    path ('porechop/*')
    path ('filtlong/*')
    path ('pycoqc/*')
    path ('kraken2_short/*')
    path ('kraken2_long/*')
    path ('quast/*')
    path ('busco/*')
    path ('prokka/*')
    path ('bakta/*')
    path ('extra/*')

    output:
    path "*multiqc_report.html"         , emit: report
    path "*_data"                       , emit: data
    path "*_assembly_metrics.csv"       , optional:true, emit: csv_assembly
    path "*_plots"                      , optional:true, emit: plots
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    """
    custom_config_args="$custom_config"

    ## Run MultiQC once to parse tool logs
    multiqc -f $args \$custom_config_args .

    ## Collect additional files to be included in the report
    if [ -d extra/ ] && compgen -G "extra/*" > /dev/null; then
        cp extra/* multiqc_data/
    fi

    ## Create the custom assembly table only when KmerFinder was run.
    if [ -s multiqc_data/multiqc_kmerfinder.yaml ]; then
        multiqc_to_custom_csv.py --assembly_type $params.assembly_type
        printf "%s\n" "exclude_modules:" "  - general_stats" > multiqc_kmerfinder_config.yaml
        custom_config_args="\$custom_config_args --config multiqc_kmerfinder_config.yaml"
    fi

    ## Run multiqc a second time
    multiqc -f $args \$custom_config_args .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html

    mkdir -p multiqc_data
    touch multiqc_data/multiqc_data.json
    touch multiqc_data/multiqc_sources.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
