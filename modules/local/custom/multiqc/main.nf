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
    path ('nanoplot/*')
    path ('porechop/*')
    path ('filtlong/*')
    path ('pycoqc/*')
    path ('kraken2_short/*')
    path ('kraken2_long/*')
    path ('quast/*')
    path ('prokka/*')
    path ('bakta/*')
    path ('extra/*')

    output:
    path "*multiqc_report.html"         , emit: report
    path "*_data"                       , emit: data
    path "*_assembly_metrics_mqc.csv"   , optional:true, emit: csv_assembly
    path "*_plots"                      , optional:true, emit: plots
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $args $custom_config .

    ## Collect additional files to be included in the report
    if [ -d extra/ ]; then
        cp extra/* multiqc_data/
    fi

    ## Create multiqc custom data
    multiqc_to_custom_csv.py --assembly_type $params.assembly_type

    ## Avoid the custom Multiqc table when the kmerfinder process is not invoked.
    if grep ">skip_kmerfinder<" workflow_summary_mqc.yaml; then
        rm *_assembly_metrics_mqc.csv
    fi

    ## Run multiqc a second time
    multiqc -f $args $custom_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
