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
      --input                       The design file used for running the pipeline in TSV format.

    Pipeline arguments:
      --assembler                   Default: "Unicycler", Available: "Canu", "Miniasm", "Unicycler". Short & Hybrid assembly always runs "Unicycler".
      --assembly_type               Default: "Short", Available: "Short", "Long", "Hybrid".
      --kraken2db                   Path to Kraken2 Database directory
      --prokka_args                 Advanced: Extra arguments to Prokka (quote and add leading space)
      --unicycler_args              Advanced: Extra arguments to Unicycler (quote and add leading space)
      --canu_args                   Advanced: Extra arguments for Canu assembly (quote and add leading space)

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
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

if(! params.skip_kraken2){
    if(params.kraken2db){
      kraken2db = file(params.kraken2db)
    } else {
      exit 1, "Missing Kraken2 DB arg"
    }
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

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
summary['Extra Canu arguments'] = params.canu_args
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

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-bacass-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/bacass Workflow Summary'
    section_href: 'https://github.com/nf-core/bacass'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


//Check compatible parameters
if(("${params.assembler}" == 'canu' || "${params.assembler}" == 'miniasm') && ("${params.assembly_type}" == 'short' || "${params.assembly_type}" == 'hybrid')){
    exit 1, "Canu and Miniasm can only be used for long read assembly and neither for Hybrid nor Shortread assembly!"
}


/* Trim and combine short read read-pairs per sample. Similar to nf-core vipr
 */
process trim_and_combine {
    label 'medium'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimming/shortreads/", mode: params.publish_dir_mode

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
    label 'medium'
    publishDir "${params.outdir}/${sample_id}/trimming/longreads/", mode: params.publish_dir_mode

    when: params.assembly_type == 'hybrid' || params.assembly_type == 'long'

    input:
    set sample_id, file(lr) from ch_for_long_trim

    output:
    set sample_id, file('trimmed.fastq') into (ch_long_trimmed_unicycler, ch_long_trimmed_canu, ch_long_trimmed_miniasm, ch_long_trimmed_consensus, ch_long_trimmed_nanopolish, ch_long_trimmed_kraken, ch_long_trimmed_medaka)
    file ("porechop.version.txt") into ch_porechop_version

    when: !('short' in params.assembly_type)

    script:
    """
    porechop -i "${lr}" -t "${task.cpus}" -o trimmed.fastq
    porechop --version > porechop.version.txt
    """
}

/*
 * STEP 1 - FastQC FOR SHORT READS
*/
process fastqc {
    label 'small'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/FastQC", mode: params.publish_dir_mode

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
    label 'medium'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/QC_longreads/NanoPlot", mode: params.publish_dir_mode

    when: (params.assembly_type != 'short')

    input:
    set sample_id, file(lr) from ch_for_nanoplot 

    output:
    file '*.png'
    file '*.html'
    file '*.txt'
    file 'nanoplot.version.txt' into ch_nanoplot_version

    script:
    """
    NanoPlot -t "${task.cpus}" --title "${sample_id}" -c darkblue --fastq ${lr}
    NanoPlot --version | sed -e "s/NanoPlot //g" > nanoplot.version.txt
    """
}


/** Quality check for nanopore Fast5 files
*/

process pycoqc{
    label 'medium'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/QC_longreads/PycoQC", mode: params.publish_dir_mode

    when: (params.assembly_type == 'hybrid' || params.assembly_type == 'long') && !params.skip_pycoqc && fast5

    input:
    set sample_id, file(lr), file(fast5) from ch_for_pycoqc.join(ch_fast5_for_pycoqc)

    output:
    set sample_id, file('sequencing_summary.txt') into ch_summary_index_for_nanopolish
    file("pycoQC_${sample_id}*")
    file("pycoQC.version.txt") into ch_pycoqc_version

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
    pycoQC --version | sed -e "s/pycoQC v//g" > pycoQC.version.txt
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
    label 'large'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/unicycler", mode: params.publish_dir_mode

    when: params.assembler == 'unicycler'

    input:
    set sample_id, file(fq1), file(fq2), file(lrfastq) from ch_short_long_joint_unicycler 

    output:
    set sample_id, file("${sample_id}_assembly.fasta") into (quast_ch, prokka_ch, dfast_ch)
    file("${sample_id}_assembly.fasta") into (ch_assembly_nanopolish_unicycler,ch_assembly_medaka_unicycler)
    file("${sample_id}_assembly.gfa")
    file("${sample_id}_unicycler.log")
    file("unicycler.version.txt") into ch_unicycler_version
    
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
    unicycler --version | sed -e "s/Unicycler v//g" > unicycler.version.txt
    """
}

process miniasm_assembly {
    label 'large'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/miniasm", mode: params.publish_dir_mode, pattern: 'assembly.fasta'
    
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
    label 'large'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/miniasm/consensus", mode: params.publish_dir_mode, pattern: 'assembly_consensus.fasta'

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
    label 'large'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/canu", mode: params.publish_dir_mode, pattern: 'assembly.fasta'

    input:
    set sample_id, file(lrfastq), val(genomeSize) from ch_long_trimmed_canu.join(ch_genomeSize_forCanu)
    
    output:
    file 'assembly.fasta' into (assembly_from_canu_for_nanopolish, assembly_from_canu_for_medaka)
    file 'canu.version.txt' into ch_canu_version

    when: params.assembler == 'canu'

    script:
    """
    canu -p assembly -d canu_out \
        genomeSize="${genomeSize}" -nanopore "${lrfastq}" \
        maxThreads="${task.cpus}" merylMemory="${task.memory.toGiga()}G" \
        merylThreads="${task.cpus}" hapThreads="${task.cpus}" batMemory="${task.memory.toGiga()}G" \
        redMemory="${task.memory.toGiga()}G" redThreads="${task.cpus}" \
        oeaMemory="${task.memory.toGiga()}G" oeaThreads="${task.cpus}" \
        corMemory="${task.memory.toGiga()}G" corThreads="${task.cpus}" ${params.canu_args}
    mv canu_out/assembly.contigs.fasta assembly.fasta
    canu --version | sed -e "s/Canu //g" > canu.version.txt
    """
}

/* kraken classification: QC for sample purity, only short end reads for now
 */
process kraken2 {
    label 'large'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/kraken", mode: params.publish_dir_mode

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
    publishDir "${params.outdir}/${sample_id}/kraken_long", mode: params.publish_dir_mode

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
  label 'small'
  tag {"$sample_id"}
  publishDir "${params.outdir}/${sample_id}/QUAST", mode: params.publish_dir_mode
  
  input:
  set sample_id, file(fasta) from quast_ch
  
  output:
  // multiqc only detects a file called report.tsv. to avoid
  // name clash with other samples we need a directory named by sample
  file("${sample_id}_assembly_QC/")
  file("${sample_id}_assembly_QC/${sample_id}_report.tsv") into quast_logs_ch
  file("quast.version.txt") into ch_quast_version

  script:
  """
  quast -t ${task.cpus} -o ${sample_id}_assembly_QC ${fasta}
  quast --version | sed -e "s/QUAST v//g" > quast.version.txt
  mv ${sample_id}_assembly_QC/report.tsv ${sample_id}_assembly_QC/${sample_id}_report.tsv
  """
}

/*
 * Annotation with prokka
 */
process prokka {
   label 'large'
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: params.publish_dir_mode
   
   input:
   set sample_id, file(fasta) from prokka_ch

   output:
   file("${sample_id}_annotation/")
   file("prokka.version.txt") into ch_prokka_version

   when: !params.skip_annotation && params.annotation_tool == 'prokka'

   script:
   """
   prokka --cpus ${task.cpus} --prefix "${sample_id}" --outdir ${sample_id}_annotation ${params.prokka_args} ${fasta}
   prokka --version | sed -e "s/prokka //g" > prokka.version.txt
   """
}

process dfast {
   label 'medium_extramem' 
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: params.publish_dir_mode

   input:
   set sample_id, file(fasta) from dfast_ch
   file (config) from Channel.value(params.dfast_config ? file(params.dfast_config) : "")

   output:
   file("RESULT*")
   file("dfast.version.txt") into ch_dfast_version

   when: !params.skip_annotation && params.annotation_tool == 'dfast'

   script:
   """
   dfast --genome ${fasta} --config $config
   dfast --version | sed -e "s/DFAST ver. //g" > dfast.version.txt
   """
}


//Polishes assembly using FAST5 files
process nanopolish {
    tag "$assembly"
    label 'large'

    publishDir "${params.outdir}/${sample_id}/nanopolish/", mode: params.publish_dir_mode, pattern: 'polished_genome.fa'

    input:
    file(assembly) from ch_assembly_consensus_for_nanopolish.mix(ch_assembly_nanopolish_unicycler,assembly_from_canu_for_nanopolish) //Should take either miniasm, canu, or unicycler consensus sequence (!)
    set sample_id, file(lrfastq), file(fast5) from ch_long_trimmed_nanopolish.join(ch_fast5_for_nanopolish)

    output:
    file 'polished_genome.fa'
    file 'nanopolish.version.txt' into ch_nanopolish_version
    file 'samtools.version.txt' into ch_samtools_version

    when: !params.skip_polish && params.assembly_type == 'long' && params.polish_method != 'medaka'

    script:
    """
    nanopolish index -d "${fast5}" "${lrfastq}"
    minimap2 -ax map-ont -t ${task.cpus} "${assembly}" "${lrfastq}"| \
    samtools sort -o reads.sorted.bam -T reads.tmp -
    samtools index reads.sorted.bam
    nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r "${lrfastq}" -b reads.sorted.bam -g "${assembly}" -t "${task.cpus}" --min-candidate-frequency 0.1
    nanopolish vcf2fasta -g "${assembly}" polished.*.vcf > polished_genome.fa

    #Versions
    nanopolish --version | sed -e "s/nanopolish version //g" | head -n 1 > nanopolish.version.txt
    samtools --version | sed -e "s/samtools //g" | head -n 1 > samtools.version.txt
    """
}

//Polishes assembly
process medaka {
    tag "$assembly"
    label 'large'

    publishDir "${params.outdir}/${sample_id}/medaka/", mode: params.publish_dir_mode, pattern: 'polished_genome.fa'

    input:
    file(assembly) from ch_assembly_consensus_for_medaka.mix(ch_assembly_medaka_unicycler,assembly_from_canu_for_medaka) //Should take either miniasm, canu, or unicycler consensus sequence (!)
    set sample_id, file(lrfastq) from ch_long_trimmed_medaka

    output:
    file 'polished_genome.fa'
    file 'medaka.version.txt' into ch_medaka_version

    when: !params.skip_polish && params.assembly_type == 'long' && params.polish_method == 'medaka'

    script:
    """
    medaka_consensus -i ${lrfastq} -d ${assembly} -o "polished_genome.fa" -t ${task.cpus}
    medaka --version | sed -e "s/medaka //g" > medaka.version.txt
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    input:
    path quast_version from ch_quast_version.first().ifEmpty([])
    path porechop_version from ch_porechop_version.first().ifEmpty([])
    path pycoqc_version from ch_pycoqc_version.first().ifEmpty([])
    path unicycler_version from ch_unicycler_version.first().ifEmpty([])
    path canu_version from ch_canu_version.first().ifEmpty([])
    path prokka_version from ch_prokka_version.first().ifEmpty([])
    path dfast_version from ch_dfast_version.first().ifEmpty([])
    path nanopolish_version from ch_nanopolish_version.first().ifEmpty([])
    path samtools_version from ch_samtools_version.first().ifEmpty([])
    path nanoplot_version from ch_nanoplot_version.first().ifEmpty([])
    path medaka_version from ch_medaka_version.first().ifEmpty([])

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    #All in main container
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    fastqc --version | sed -e "s/FastQC v//g" > fastqc.version.txt
    
    #Inside main container
    miniasm -V > miniasm.version.txt
    minimap2 --version &> minimap2.version.txt
    racon --version | sed -e "s/v//g" > racon.version.txt
    skewer --version | sed -e "s/skewer version://g" | sed -e 's/\\s//g' | head -n 1  > skewer.version.txt
    kraken2 --version | sed -e "s/Kraken version //g" | head -n 1 > kraken2.version.txt
    multiqc --version | sed -e "s/multiqc, version//g" > multiqc.version.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * STEP - MultiQC
 */

process multiqc {
    label 'small'
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:

    path (multiqc_config) from ch_multiqc_config
    path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    
    //file prokka_logs from prokka_logs_ch.collect().ifEmpty([])
    file ('quast_logs/*') from quast_logs_ch.collect().ifEmpty([])
    // NOTE unicycler and kraken not supported
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/bacass] Successful: $workflow.runName"
    if (!workflow.success) {
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
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/bacass] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/bacass] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
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
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/bacass] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/bacass] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/bacass]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/bacass]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/bacass v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
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
