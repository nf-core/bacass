process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path software_versions
   //path workflow_summary
    path multiqc_logo
    path ('fastqc/*')
    path ('fastp/*')
    path ('nanoplot/*')
    path ('porechop/*')
    path ('pycoqc/*')
    path ('kraken2_short/*')
    path ('kraken2_long/*')
    path ('quast_unicycler/*')
    path ('prokka/*')
    path ('bakta/*')
    path ('extra/*')

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path "*_plots"                  , optional:true, emit: plots
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $args $custom_config .

    ## Collect extra fields to be included in the report
    cp extra/* multiqc_data/

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py

    ## Run multiqc a second time
    multiqc -f $args -e general_stats $custom_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
