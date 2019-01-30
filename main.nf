#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/bacass
========================================================================================
 nf-core/bacass Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/bacass
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/bacass v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/bacass -params-file params.yaml -profile docker
    or
    nextflow run nf-core/bacass --reads '*_R{1,2}.fastq.gz' --kraken2db 'path-to-kraken2db' -profile docker

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Other options:
      -params-file                  A parameters file, listing parameters defined here, including paired FastQ input
				    if not defined through `--reads` (see below). See docs for more info
      --reads                       Path to paired-end input reads (must be surrounded with quotes; see also docs)
                                    Can be replaced with samples dictionary in params-file (see above)
      --skip_kraken2                Don't run Kraken2 for classification
      --kraken2db                   Path to Kraken2 Database directory
      --unicycler_args              Advanced: Extra arguments to Unicycler (quote and add leading space)
      --prokka_args                 Advanced: Extra arguments to Prokka (quote and add leading space)
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// see https://ccb.jhu.edu/software/kraken2/index.shtml#downloads
if(!params.skip_kraken2) {
    if (!params.kraken2db)
        exit 1, "Missing Kraken2 DB arg"
    kraken2db = file(params.kraken2db)
    if (!kraken2db.exists())
        exit 1, "Missing Kraken2 DB: ${kraken2db}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/* GetReadPair and GetReadUnitKeys support read pair input via yaml files
 * This allows for merging of libraries for samples that were sequenced twice
 * and is generally better for reproducibility and keeping of records (all
 * config done via files)
 */
def GetReadPair = { sk, rk ->
    // FIXME if files don't exist, their path might be relative to the input yaml
    // see https://gist.github.com/ysb33r/5804364
    tuple(file(params.samples[sk].readunits[rk]['fq1']),
          file(params.samples[sk].readunits[rk]['fq2']))
}

def GetReadUnitKeys = { sk ->
    params.samples[sk].readunits.keySet()
}

if (params.reads) {
    Channel.fromFilePairs( params.reads )// flat: true
	.set { fastq_ch }

} else if (params.readPaths) {
   Channel.from( params.readPaths )
	.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
        .set { fastq_ch }

} else {
    sample_keys = params.samples? params.samples.keySet() : []
    Channel.from( sample_keys )
        .map { sk -> tuple(sk, GetReadUnitKeys(sk).collect{GetReadPair(sk, it)}.flatten()) }
        .set { fastq_ch }
}
//println "List of samples: " +  sample_keys.join(", ")
//fastq_ch.subscribe { println "$it" }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/bacass v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/bacass'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name'] = custom_runName ?: workflow.runName
//summary['Sample keys'] = sample_keys
summary['Skip Kraken2'] = params.skip_kraken2
summary['Kraken2 DB'] = params.kraken2db
summary['Extra Unicycler arguments'] = params.unicycler_args
summary['Extra Prokka arguments'] = params.prokka_args
summary['Max Memory'] = params.max_memory
summary['Max CPUs'] = params.max_cpus
summary['Max Time'] = params.max_time
summary['Output dir'] = params.outdir
summary['Working dir'] = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-bacass-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/bacass Workflow Summary'
    section_href: 'https://github.com/nf-core/bacass'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    prokka -v 2> v_prokka.txt
    quast -v > v_quast.txt
    skewer -v > v_skewer.txt
    kraken2 -v > v_kraken2.txt
    Bandage -v > v_bandage.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



/* Trim and combine read-pairs per sample. Similar to nf-core vipr
 */
process trim_and_combine {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/${sample_id}_reads/", mode: 'copy'

    input:
    set sample_id, file(reads) from fastq_ch

    output:
    set sample_id, file("${sample_id}_trm-cmb.R1.fastq.gz"), file("${sample_id}_trm-cmb.R2.fastq.gz") \
	into kraken2_ch, unicycler_ch, fastqc_ch
    // not keeping logs for multiqc input. for that to be useful we would need to concat first and then run skewer
    
    script:
    """
    # loop over readunits in pairs per sample
    pairno=0
    echo ${reads.join(" ")} | xargs -n2 | while read fq1 fq2; do
	skewer --quiet -t ${task.cpus} -m pe -q 3 -n -z \$fq1 \$fq2;
    done
    cat \$(ls *trimmed-pair1.fastq.gz | sort) >> ${sample_id}_trm-cmb.R1.fastq.gz
    cat \$(ls *trimmed-pair2.fastq.gz | sort) >> ${sample_id}_trm-cmb.R2.fastq.gz
    """
}


/* fastqc
 */
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/${sample_id}_reads", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from fastqc_ch

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -t {task.cpus} -q ${fq1} ${fq2}
    """
}


/* unicycler
 */
process unicycler {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from unicycler_ch

    output:
    set sample_id, file("${sample_id}_assembly.fasta") into quast_ch, prokka_ch
    set sample_id, file("${sample_id}_assembly.gfa") into bandage_ch
    file("${sample_id}_assembly.gfa")
    file("${sample_id}_assembly.png")
    file("${sample_id}_unicycler.log")
    
    script:
    """
    unicycler -1 $fq1 -2 $fq2 --threads ${task.cpus} ${params.unicycler_args} --keep 0 -o .
    mv unicycler.log ${sample_id}_unicycler.log
    # rename so that quast can use the name 
    mv assembly.gfa ${sample_id}_assembly.gfa
    mv assembly.fasta ${sample_id}_assembly.fasta
    Bandage image ${sample_id}_assembly.gfa ${sample_id}_assembly.png
    """
}


/* waste to have a separate process because it's superfast
process bandage {
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

   input:
   set sample_id, gfa from bandage_ch

   output:
   file("${sample_id}_assembly_bandage.png")

   script:
   """
   Bandage image ${gfa} ${sample_id}_assembly_bandage.png
   """
}
*/


/* kraken classification: QC for sample purity
 */
if(!params.skip_kraken2) {
    process kraken2 {
        tag "$sample_id"
        publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

        input:
        set sample_id, file(fq1), file(fq2) from kraken2_ch

        output:
        file("${sample_id}_kraken2.report")

        script:
	"""
        # stdout reports per read which is not needed. kraken.report can be used with pavian
        # braken would be nice but requires readlength and correspondingly build db
	kraken2 --threads ${task.cpus} --paired --db ${kraken2db} \
		--report ${sample_id}_kraken2.report ${fq1} ${fq2} | gzip > kraken2.out.gz
	"""
    }
}


/* assembly qc with quast
 */
process quast {
  tag { "quast for each $sample_id" }
  publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
  
  input:
  set sample_id, fasta from quast_ch
  
  output:
  // multiqc only detects a file called report.tsv. to avoid
  // name clash with other samples we need a directory named by sample
  file("${sample_id}_assembly_QC/") into quast_logs_ch

  script:
  """
  quast.py -t ${task.cpus} -o ${sample_id}_assembly_QC ${fasta} 
  """
}


/* annotation with prokka
 */
process prokka {
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

   input:
   set sample_id, fasta from prokka_ch

   output:
   file("${sample_id}_annotation/")
   // multiqc prokka module is just a stub using txt. see https://github.com/ewels/MultiQC/issues/587
   // also, this only makes sense if we could set genus/species/strain. otherwise all samples
   // are the same
   // file("${sample_id}_annotation/*txt") into prokka_logs_ch

   script:
   """
   prokka --cpus ${task.cpus} --prefix "${sample_id}" --outdir ${sample_id}_annotation ${params.prokka_args} ${fasta}
   """
}


/* Multiqc
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    //file prokka_logs from prokka_logs_ch.collect().ifEmpty([])
    file quast_logs from quast_logs_ch.collect().ifEmpty([])
    // NOTE unicycler and kraken not supported
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


/* Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


// to enable poor man's syntax checking
def testSyntaxOnly(){}


/*
 * Completion e-mail notification
 */
 workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/bacass] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/bacass] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/bacass] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/bacass] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/bacass] Pipeline Complete"

}

