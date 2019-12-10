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
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/bacass --input input.csv --kraken2db 'path-to-kraken2db' -profile docker

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --input                       The design file used for running the pipeline.

    Pipeline arguments:
        --assembler                   Default: "Unicycler", Available: "Canu", "Miniasm", "Unicycler". Short & Hybrid assembly always runs "Unicycler".
        --assembly_type               Default: "Short", Available: "Short", "Long", "Hybrid".
        --kraken2db                   Path to Kraken2 Database directory
        --prokka_args                 Advanced: Extra arguments to Prokka (quote and add leading space)
        --unicycler_args              Advanced: Extra arguments to Unicycler (quote and add leading space)
  
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
   Skipping options:
      --skip_annotation             Skips the annotation with Prokka
      --skip_kraken2                Skips the read classification with Kraken2
      --skip_polish                 Skips polishing long-reads with Nanopolish or Medaka
      --skip_pycoqc                 Skips long-read raw signal QC

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


if(!params.skip_kraken2){
    if(params.kraken2db){
      kraken2db = file(params.kraken2db)
    } else {
      exit 1, "Missing Kraken2 DB arg"
    }
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
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")


//Check whether we have a design file as input set
if(!params.input){
    exit 1, "Missing Design File - please see documentation how to create one."
} else {
    //Design file looks like this
    // ID R1 R2 Long-ReadFastQ Fast5Path GenomeSize
    // ID is required, everything else (can!) be optional and causes some pipeline components to turn off!
    // Tapping the parsed input design to multiple channels to get some data to specific downstream processes that don't need full information!
    Channel
    .fromPath(params.input)
    .splitCsv(header: true, sep:'\t')
    .map { col -> 
           def id = "${col.ID}" 
           def r1 = returnFile("${col.R1}")
           def r2 = returnFile("${col.R2}")
           def lr = returnFile("${col.LongFastQ}")
           def f5 = returnFile("${col.Fast5}")
           def genome_size = "${col.GenomeSize}"
           tuple(id,r1,r2,lr,f5,genome_size)
    }
    .dump(tag: "input")
    .tap {ch_all_data; ch_all_data_for_fast5; ch_all_data_for_genomesize}
    .map { id,r1,r2,lr,f5,gs -> 
    tuple(id,r1,r2) 
    }
    .filter{ id,r1,r2 -> 
    r1 != 'NA' && r2 != 'NA'}
    //Filter to get rid of R1/R2 that are NA
    .into {ch_for_short_trim; ch_for_fastqc}
    //Dump long read info to different channel! 
    ch_all_data
    .map { id, r1, r2, lr, f5, genomeSize -> 
            tuple(id, file(lr))
    }
    .dump(tag: 'longinput')
    .into {ch_for_long_trim; ch_for_nanoplot; ch_for_pycoqc; ch_for_nanopolish; ch_for_long_fastq}

    //Dump fast5 to separate channel
    ch_all_data_for_fast5
    .map { id, r1, r2, lr, f5, genomeSize -> 
            tuple(id, f5)
    }
    .filter {id, fast5 -> 
        fast5 != 'NA'
    }
    .into {ch_fast5_for_pycoqc; ch_fast5_for_nanopolish}

    //Dump genomeSize to separate channel, too
    ch_all_data_for_genomesize
    .map { id, r1, r2, lr, f5, genomeSize -> 
    tuple(id,genomeSize)
    }
    .filter{id, genomeSize -> 
      genomeSize != 'NA'
    }
    .set {ch_genomeSize_forCanu}
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Pipeline Name'] = 'nf-core/bacass'
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Assembler Method'] = params.assembler
summary['Assembly Type'] = params.assembly_type
if (params.kraken2db) summary['Kraken2 DB'] = params.kraken2db 
summary['Extra Prokka arguments'] = params.prokka_args
summary['Extra Unicycler arguments'] = params.unicycler_args
if (params.skip_annotation) summary['Skip Annotation'] = params.skip_annotation
if (params.skip_kraken2) summary['Skip Kraken2'] = params.skip_kraken2
if (params.skip_polish) summary['Skip Polish'] = params.skip_polish
if (!params.skip_polish) summary['Polish Method'] = params.polish_method
if (params.skip_pycoqc) summary['Skip PycoQC'] = params.skip_pycoqc
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Launch dir']       = workflow.launchDir
summary['Output dir'] = params.outdir
summary['Working dir'] = workflow.workDir
summary['Script dir'] = workflow.projectDir
summary['User'] = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
summary['Config Profile'] = workflow.profile

if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"


// Check the hostnames against configured profiles
checkHostname()

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

//Check compatible parameters
if(("${params.assembler}" == 'canu' || "${params.assembler}" == 'miniasm') && ("${params.assembly_type}" == 'short' || "${params.assembly_type}" == 'hybrid')){
    exit 1, "Canu and Miniasm can only be used for long read assembly and neither for Hybrid nor Shortread assembly!"
}


/* Trim and combine short read read-pairs per sample. Similar to nf-core vipr
 */
process trim_and_combine {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimming/shortreads/", mode: 'copy'

    label 'medium'

    input:
    set sample_id, file(r1), file(r2) from ch_for_short_trim

    output:
    set sample_id, file("${sample_id}_trm-cmb.R1.fastq.gz"), file("${sample_id}_trm-cmb.R2.fastq.gz") into (ch_short_for_kraken2, ch_short_for_unicycler, ch_short_for_fastqc)
    // not keeping logs for multiqc input. for that to be useful we would need to concat first and then run skewer
    
    script:
    """
    # loop over readunits in pairs per sample
    pairno=0
    echo "${r1} ${r2}" | xargs -n2 | while read fq1 fq2; do
	skewer --quiet -t ${task.cpus} -m pe -q 3 -n -z \$fq1 \$fq2;
    done
    cat \$(ls *trimmed-pair1.fastq.gz | sort) >> ${sample_id}_trm-cmb.R1.fastq.gz
    cat \$(ls *trimmed-pair2.fastq.gz | sort) >> ${sample_id}_trm-cmb.R2.fastq.gz
    """
}


//AdapterTrimming for ONT reads
process adapter_trimming {
    
    publishDir "${params.outdir}/${sample_id}/trimming/longreads/", mode: 'copy'

    when: params.assembly_type == 'hybrid' || params.assembly_type == 'long'

    label 'medium'

    input:
	set sample_id, file(lr) from ch_for_long_trim

    output:
    set sample_id, file('trimmed.fastq') into (ch_long_trimmed_unicycler, ch_long_trimmed_canu, ch_long_trimmed_miniasm, ch_long_trimmed_consensus, ch_long_trimmed_nanopolish, ch_long_trimmed_kraken, ch_long_trimmed_medaka)
    file ("v_porechop.txt") into ch_porechop_version

    when: !('short' in params.assembly_type)


	script:
    """
    porechop -i "${lr}" -t "${task.cpus}" -o trimmed.fastq
    porechop --version > v_porechop.txt
    """
}

/*
 * STEP 1 - FastQC FOR SHORT READS
*/
process fastqc {
    label 'small'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/FastQC", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from ch_short_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -q ${fq1} ${fq2}
    """
}

/*
 * Quality check for nanopore reads and Quality/Length Plots
 */
process nanoplot {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/QC_longreads/NanoPlot", mode: 'copy'

    when: (params.assembly_type != 'short')

    input:
    set sample_id, file(lr) from ch_for_nanoplot 

    output:
    file '*.png'
    file '*.html'
    file '*.txt'

    script:
    """
    NanoPlot -t "${task.cpus}" --title "${sample_id}" -c darkblue --fastq ${lr}
    """
}


/** Quality check for nanopore Fast5 files
*/

process pycoqc{
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/QC_longreads/PycoQC", mode: 'copy'

    when: (params.assembly_type == 'hybrid' || params.assembly_type == 'long') && !params.skip_pycoqc && fast5

    input:
    set sample_id, file(lr), file(fast5) from ch_for_pycoqc.join(ch_fast5_for_pycoqc)

    output:
    set sample_id, file('sequencing_summary.txt') into ch_summary_index_for_nanopolish
    file("pycoQC_${sample_id}*")

    script:
    //Find out whether the sequencing_summary already exists
    if(file("${fast5}/sequencing_summary.txt").exists()){
        run_summary = ''
        prefix = "${fast5}/"
    } else {
        run_summary =  "Fast5_to_seq_summary -f $fast5 -t ${task.cpus} -s './sequencing_summary.txt' --verbose_level 2"
        prefix = ''
    }
    //Barcodes available? 
    barcode_me = file("${fast5}/barcoding_sequencing.txt").exists() ? "-b ${fast5}/barcoding_sequencing.txt" : ''
    """
    $run_summary
    pycoQC -f "${prefix}sequencing_summary.txt" $barcode_me -o pycoQC_${sample_id}.html -j pycoQC_${sample_id}.json
    """
}

/* Join channels for unicycler, as trimming the files happens in two separate processes for paralellization of individual steps. As samples have the same sampleID, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
*/ 
if(params.assembly_type == 'hybrid'){
    ch_short_for_unicycler
        .join(ch_long_trimmed_unicycler)
        .dump(tag: 'unicycler')
        .set {ch_short_long_joint_unicycler}
} else if(params.assembly_type == 'short'){
    ch_short_for_unicycler
        .map{id,R1,R2 -> 
        tuple(id,R1,R2,'NA')}
        .dump(tag: 'unicycler')
        .set {ch_short_long_joint_unicycler}
} else if(params.assembly_type == 'long'){
    ch_long_trimmed_unicycler
        .map{id,lr -> 
        tuple(id,'NA','NA',lr)}
        .dump(tag: 'unicycler')
        .set {ch_short_long_joint_unicycler}
}

/* unicycler (short, long or hybrid mode!)
 */
process unicycler {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/unicycler", mode: 'copy'

    when: params.assembler == 'unicycler'

    input:
    set sample_id, file(fq1), file(fq2), file(lrfastq) from ch_short_long_joint_unicycler 

    output:
    set sample_id, file("${sample_id}_assembly.fasta") into (quast_ch, prokka_ch, dfast_ch)
    set sample_id, file("${sample_id}_assembly.gfa") into bandage_ch
    file("${sample_id}_assembly.fasta") into (ch_assembly_nanopolish_unicycler,ch_assembly_medaka_unicycler)
    file("${sample_id}_assembly.gfa")
    file("${sample_id}_assembly.png")
    file("${sample_id}_unicycler.log")
    
    script:
    if(params.assembly_type == 'long'){
        data_param = "-l $lrfastq"
    } else if (params.assembly_type == 'short'){
        data_param = "-1 $fq1 -2 $fq2"
    } else if (params.assembly_type == 'hybrid'){
        data_param = "-1 $fq1 -2 $fq2 -l $lrfastq"
    }

    """
    unicycler $data_param --threads ${task.cpus} ${params.unicycler_args} --keep 0 -o .
    mv unicycler.log ${sample_id}_unicycler.log
    # rename so that quast can use the name 
    mv assembly.gfa ${sample_id}_assembly.gfa
    mv assembly.fasta ${sample_id}_assembly.fasta
    Bandage image ${sample_id}_assembly.gfa ${sample_id}_assembly.png
    """
}

process miniasm_assembly {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/miniasm", mode: 'copy', pattern: 'assembly.fasta'
    
    label 'large'


    input:
    set sample_id, file(lrfastq) from ch_long_trimmed_miniasm

    output:
    file 'assembly.fasta' into ch_assembly_from_miniasm

    when: params.assembler == 'miniasm'

    script:
    """
    minimap2 -x ava-ont -t "${task.cpus}" "${lrfastq}" "${lrfastq}" > "${lrfastq}.paf"
    miniasm -f "${lrfastq}" "${lrfastq}.paf" > "${lrfastq}.gfa"
    awk '/^S/{print ">"\$2"\\n"\$3}' "${lrfastq}.gfa" | fold > assembly.fasta
    """
}

//Run consensus for miniasm, the others don't need it.
process consensus {
    tag "$sample_id"
	publishDir "${params.outdir}/${sample_id}/miniasm/consensus", mode: 'copy', pattern: 'assembly_consensus.fasta'
    label 'large'

    input:
    set sample_id, file(lrfastq) from ch_long_trimmed_consensus
    file(assembly) from ch_assembly_from_miniasm

    output:
    file 'assembly_consensus.fasta' into (ch_assembly_consensus_for_nanopolish, ch_assembly_consensus_for_medaka)

	script:
    """
    minimap2 -x map-ont -t "${task.cpus}" "${assembly}" "${lrfastq}" > assembly.paf
    racon -t "${task.cpus}" "${lrfastq}" assembly.paf "${assembly}" > assembly_consensus.fasta
    """
}

process canu_assembly {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/canu", mode: 'copy', pattern: 'assembly.fasta'

    label 'large'

    input:
    set sample_id, file(lrfastq), val(genomeSize) from ch_long_trimmed_canu.join(ch_genomeSize_forCanu)
    
    output:
    file 'assembly.fasta' into (assembly_from_canu_for_nanopolish, assembly_from_canu_for_medaka)

    when: params.assembler == 'canu'

    script:
    """
    canu -p assembly -d canu_out \
        genomeSize="${genomeSize}" -nanopore-raw "${lrfastq}" \
        maxThreads="${task.cpus}" merylMemory="${task.memory.toGiga()}G" \
        merylThreads="${task.cpus}" hapThreads="${task.cpus}" batMemory="${task.memory.toGiga()}G" \
        corMemory="${task.memory.toGiga()}G" corThreads="${task.cpus}" ${params.canu_args}
    mv canu_out/assembly.contigs.fasta assembly.fasta
    """
}

/* kraken classification: QC for sample purity, only short end reads for now
 */
process kraken2 {
    label 'large'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/kraken", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from ch_short_for_kraken2

    output:
    file("${sample_id}_kraken2.report")

    when: !params.skip_kraken2

    script:
	"""
    # stdout reports per read which is not needed. kraken.report can be used with pavian
    # braken would be nice but requires readlength and correspondingly build db
	kraken2 --threads ${task.cpus} --paired --db ${kraken2db} \
		--report ${sample_id}_kraken2.report ${fq1} ${fq2} | gzip > kraken2.out.gz
	"""
}

/* kraken classification: QC for sample purity, only short end reads for now
 */
process kraken2_long {
    label 'large'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/kraken_long", mode: 'copy'

    input:
    set sample_id, file(lr) from ch_long_trimmed_kraken

    output:
    file("${sample_id}_kraken2.report")

    when: !params.skip_kraken2

    script:
	"""
    # stdout reports per read which is not needed. kraken.report can be used with pavian
    # braken would be nice but requires readlength and correspondingly build db
	kraken2 --threads ${task.cpus} --db ${kraken2db} \
		--report ${sample_id}_kraken2.report ${lr} | gzip > kraken2.out.gz
	"""
}

/* assembly qc with quast
 */
process quast {
  tag {"$sample_id"}
  publishDir "${params.outdir}/${sample_id}/QUAST", mode: 'copy'
  
  input:
  set sample_id, file(fasta) from quast_ch
  
  output:
  // multiqc only detects a file called report.tsv. to avoid
  // name clash with other samples we need a directory named by sample
  file("${sample_id}_assembly_QC/")
  file("${sample_id}_assembly_QC/report.tsv") into quast_logs_ch
  file("v_quast.txt") into ch_quast_version

  script:
  """
  quast -t ${task.cpus} -o ${sample_id}_assembly_QC ${fasta}
  quast -v > v_quast.txt
  """
}

/*
 * Annotation with prokka
 */
process prokka {
   label 'large'
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
   
   input:
   set sample_id, file(fasta) from prokka_ch

   output:
   file("${sample_id}_annotation/")
   // multiqc prokka module is just a stub using txt. see https://github.com/ewels/MultiQC/issues/587
   // also, this only makes sense if we could set genus/species/strain. otherwise all samples
   // are the same
   // file("${sample_id}_annotation/*txt") into prokka_logs_ch

   when: !params.skip_annotation && params.annotation_tool == 'prokka'

   script:
   """
   prokka --cpus ${task.cpus} --prefix "${sample_id}" --outdir ${sample_id}_annotation ${params.prokka_args} ${fasta}
   """
}

process dfast {
   label 'large'
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
   

   input:
   set sample_id, file(fasta) from dfast_ch

   output:
   file("${sample_id}_annotation/")
   file("v_dfast.txt") into ch_dfast_version_for_multiqc

   when: !params.skip_annotation && params.annotation_tool == 'dfast'

   script:
   """
   dfast --genome ${fasta} --config ${params.dfast_config}
   dfast &> v_dfast.txt 2>&1 || true
   """
}


//Polishes assembly using FAST5 files
process nanopolish {
    tag "$assembly"
    label 'large'

    publishDir "${params.outdir}/${sample_id}/nanopolish/", mode: 'copy', pattern: 'polished_genome.fa'

    input:
    file(assembly) from ch_assembly_consensus_for_nanopolish.mix(ch_assembly_nanopolish_unicycler,assembly_from_canu_for_nanopolish) //Should take either miniasm, canu, or unicycler consensus sequence (!)
    set sample_id, file(lrfastq), file(fast5) from ch_long_trimmed_nanopolish.join(ch_fast5_for_nanopolish)

    output:
    file 'polished_genome.fa'

    when: !params.skip_polish && params.assembly_type == 'long' && params.polish_method != 'medaka'

    script:
    """
    nanopolish index -d "${fast5}" "${lrfastq}"
    minimap2 -ax map-ont -t ${task.cpus} "${assembly}" "${lrfastq}"| \
    samtools sort -o reads.sorted.bam -T reads.tmp -
    samtools index reads.sorted.bam
    nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r "${lrfastq}" -b reads.sorted.bam -g "${assembly}" -t "${task.cpus}" --min-candidate-frequency 0.1
    nanopolish vcf2fasta -g "${assembly}" polished.*.vcf > polished_genome.fa
    """
}

//Polishes assembly
process medaka {
    tag "$assembly"
    label 'large'

    publishDir "${params.outdir}/${sample_id}/medaka/", mode: 'copy', pattern: 'polished_genome.fa'

    input:
    file(assembly) from ch_assembly_consensus_for_medaka.mix(ch_assembly_medaka_unicycler,assembly_from_canu_for_medaka) //Should take either miniasm, canu, or unicycler consensus sequence (!)
    set sample_id, file(lrfastq) from ch_long_trimmed_medaka

    output:
    file 'polished_genome.fa'

    when: !params.skip_polish && params.assembly_type == 'long' && params.polish_method == 'medaka'

    script:
    """
    medaka_consensus -i ${lrfastq} -d ${assembly} -o "polished_genome.fa" -t ${task.cpus}
    """
}




/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    input:
    file quast_version from ch_quast_version
    file porechop_version from ch_porechop_version
    file dfast_version from ch_dfast_version_for_multiqc


    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    prokka -v 2> v_prokka.txt
    skewer -v > v_skewer.txt
    kraken2 -v > v_kraken2.txt
    Bandage -v > v_bandage.txt
    nanopolish --version > v_nanopolish.txt
    miniasm -V > v_miniasm.txt
    racon --version > v_racon.txt
    samtools --version &> v_samtools.txt 2>&1 || true
    minimap2 --version &> v_minimap2.txt
    NanoPlot --version > v_nanoplot.txt
    canu --version > v_canu.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * STEP - MultiQC
 */

process multiqc {
    label 'small'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    //file prokka_logs from prokka_logs_ch.collect().ifEmpty([])
    file ('quast_logs/*') from quast_logs_ch.collect().ifEmpty([])
    // NOTE unicycler and kraken not supported
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


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
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/bacass] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/bacass] Could not attach MultiQC report to summary email"
    }

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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
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
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/bacass]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/bacass]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple} nf-core/bacass v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

// Return file if it exists, if NA is found this gets treated as a String information
static def returnFile(it) {
    if(it == 'NA') {
        return 'NA'
    } else { 
    if (!file(it).exists()) exit 1, "Warning: Missing file in CSV file: ${it}, see --help for more information"
        return file(it)
    }
}
