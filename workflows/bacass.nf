/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowBacass.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db, params.dfast_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check krakendb
if(! params.skip_kraken2){
    if(params.kraken2db){
        kraken2db = file(params.kraken2db)
    } else {
        exit 1, "Missing Kraken2 DB arg"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { PYCOQC                    } from '../modules/local/pycoqc'
include { UNICYCLER                 } from '../modules/local/unicycler'
include { NANOPOLISH                } from '../modules/local/nanopolish'
include { MEDAKA                    } from '../modules/local/medaka'
include { KRAKEN2_DB_PREPARATION    } from '../modules/local/kraken2_db_preparation'
include { KMERFINDER                } from '../modules/local/kmerfinder'
include { KMERFINDER_SUMMARY        } from '../modules/local/kmerfinder_summary'
include { DFAST                     } from '../modules/local/dfast'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { NANOPLOT                              } from '../modules/nf-core/nanoplot/main'
include { PORECHOP_PORECHOP                     } from '../modules/nf-core/porechop/porechop/main'
include { CANU                                  } from '../modules/nf-core/canu/main'
include { MINIMAP2_ALIGN                        } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CONSENSUS  } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_POLISH     } from '../modules/nf-core/minimap2/align/main'
include { MINIASM                               } from '../modules/nf-core/miniasm/main'
include { DRAGONFLYE                            } from '../modules/nf-core/dragonflye/main'
include { RACON                                 } from '../modules/nf-core/racon/main'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2            } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_LONG       } from '../modules/nf-core/kraken2/kraken2/main'
include { QUAST                                 } from '../modules/nf-core/quast/main'
include { GUNZIP                                } from '../modules/nf-core/gunzip/main'
include { PROKKA                                } from '../modules/nf-core/prokka/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_TRIM_FASTP_FASTQC               } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { BAKTA_DBDOWNLOAD_RUN                  } from '../subworkflows/local/bakta_dbdownload_run'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BACASS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    def criteria = multiMapCriteria {
        meta, fastq_1, fastq_2, long_fastq, fast5 ->
            shortreads: fastq_1     != 'NA' ? tuple(meta, [file(fastq_1), file(fastq_2)]) : null
            longreads: long_fastq   != 'NA' ? tuple(meta, file(long_fastq))         : null
            fast5: fast5            != 'NA' ? tuple(meta, fast5)              : null
    }
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    Channel
        .fromSamplesheet('input')
        .multiMap (criteria)
        .set { ch_input }
    // reconfigure channels
    ch_input
        .shortreads
        .filter{ it != null }
        .set { ch_shortreads }
    ch_input
        .longreads
        .filter{ it != null }
        .set { ch_longreads }
    ch_input
        .fast5
        .filter{ it != null }
        .set { ch_fast5 }

    //
    // SUBWORKFLOW: Short reads QC and trim adapters
    //
    FASTQ_TRIM_FASTP_FASTQC (
        ch_shortreads,
        [],
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )
    ch_fastqc_raw_multiqc   = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip
    ch_fastqc_trim_multiqc  = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip
    ch_trim_json_multiqc    = FASTQ_TRIM_FASTP_FASTQC.out.trim_json
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions.ifEmpty(null))

    //
    // MODULE: Nanoplot, quality check for nanopore reads and Quality/Length Plots
    //
    NANOPLOT (
        ch_longreads
    )
    ch_nanoplot_txt_multiqc = NANOPLOT.out.txt
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.ifEmpty(null))

    //
    // MODULE: PYCOQC, quality check for nanopore reads and Quality/Length Plots
    //
    // TODO: Couldn't be tested. No configuration test available (lack of fast5 file or params.skip_pycoqc=false).
    ch_pycoqc_multiqc = Channel.empty()
    if ( !params.skip_pycoqc ) {
        PYCOQC (
            ch_fast5.dump(tag: 'fast5')
        )
        ch_pycoqc_multiqc = PYCOQC.out.json
        ch_versions       = ch_versions.mix(PYCOQC.out.versions.ifEmpty(null))
    }

    //
    // MODULE: PORECHOP, quality check for nanopore reads and Quality/Length Plots
    //
    ch_porechop_log_multiqc = Channel.empty()
    if ( params.assembly_type == 'hybrid' || params.assembly_type == 'long' && !('short' in params.assembly_type) ) {
        PORECHOP_PORECHOP (
            ch_longreads.dump(tag: 'longreads')
        )
        ch_porechop_log_multiqc = PORECHOP_PORECHOP.out.log
        ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions.ifEmpty(null) )
    }

    //
    // Join channels for assemblers. As samples have the same meta data, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
    // Prepare channel for Kraken2
    //
    if(params.assembly_type == 'hybrid'){
        ch_for_kraken2_short    = FASTQ_TRIM_FASTP_FASTQC.out.reads
        ch_for_kraken2_long     = PORECHOP_PORECHOP.out.reads
        FASTQ_TRIM_FASTP_FASTQC.out.reads
            .dump(tag: 'fastp')
            .join(PORECHOP_PORECHOP.out.reads)
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
    } else if ( params.assembly_type == 'short' ) {
        ch_for_kraken2_short    = FASTQ_TRIM_FASTP_FASTQC.out.reads
        ch_for_kraken2_long     = Channel.empty()
        FASTQ_TRIM_FASTP_FASTQC.out.reads
            .dump(tag: 'fastp')
            .map{ meta,reads -> tuple(meta,reads,[]) }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
    } else if ( params.assembly_type == 'long' ) {
        ch_for_kraken2_short    = Channel.empty()
        ch_for_kraken2_long     = PORECHOP_PORECHOP.out.reads
        PORECHOP_PORECHOP.out.reads
            .dump(tag: 'porechop')
            .map{ meta,lr -> tuple(meta,[],lr) }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
    }
/*
    //
    // ASSEMBLY: Unicycler, Canu, Miniasm, Dragonflye
    //
    ch_assembly = Channel.empty()

    //
    // MODULE: Unicycler, genome assembly, nf-core module allows only short, long and hybrid assembly
    //
    if ( params.assembler == 'unicycler' ) {
        UNICYCLER (
            ch_for_assembly
        )
        ch_assembly = ch_assembly.mix( UNICYCLER.out.scaffolds.dump(tag: 'unicycler') )
        ch_versions = ch_versions.mix( UNICYCLER.out.versions.ifEmpty(null) )
    }


    //
    // MODULE: Canu, genome assembly, long reads
    //
    if ( params.assembler == 'canu' ) {
        CANU (
            ch_for_assembly.map { meta, reads, lr -> tuple( meta, lr ) },
            params.canu_mode,
            ch_for_assembly.map { meta, reads, lr -> meta.genome_size }
        )
        ch_assembly = ch_assembly.mix( CANU.out.assembly.dump(tag: 'canu') )
        ch_versions = ch_versions.mix(CANU.out.versions.ifEmpty(null))
    }

    //
    // MODULE: Miniasm, genome assembly, long reads
    if ( params.assembler == 'miniasm' ) {
        MINIMAP2_ALIGN (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            [],
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.ifEmpty(null))

        ch_for_assembly
            .join(MINIMAP2_ALIGN.out.paf)
            .map { meta, sr, lr, paf-> tuple(meta, lr, paf) }
            .set { ch_for_miniasm }

        MINIASM (
            ch_for_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions.ifEmpty(null))

        MINIMAP2_CONSENSUS (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            MINIASM.out.assembly.map { meta, assembly -> assembly },
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS.out.versions.ifEmpty(null))

        ch_for_assembly
            .join(MINIASM.out.assembly)
            .join(MINIMAP2_CONSENSUS.out.paf)
            .map { meta, sr, lr, assembly, paf -> tuple(meta, lr, assembly, paf) }
            .set{ ch_for_racon }

        RACON (
            ch_for_racon
        )
        ch_assembly = ch_assembly.mix( RACON.out.improved_assembly.dump(tag: 'miniasm') )
        ch_versions = ch_versions.mix( RACON.out.versions.ifEmpty(null) )
    }

    //
    // MODULE: Dragonflye, genome assembly of long reads. Moreover, it provides the option for polishing the draft genome using short reads when both short and long reads are available.
    //
    if( params.assembler == 'dragonflye' ){
        DRAGONFLYE(
            ch_for_assembly
        )
        ch_assembly = ch_assembly.mix( DRAGONFLYE.out.contigs.dump(tag: 'dragonflye') )
        ch_versions = ch_versions.mix( DRAGONFLYE.out.versions.ifEmpty(null) )
    }

    //
    // MODULE: Nanopolish, polishes assembly using FAST5 files - should take either miniasm, canu, or unicycler consensus sequence
    //
    if ( !params.skip_polish && params.assembly_type == 'long' && params.polish_method != 'medaka' ) {
        ch_for_assembly
            .join( ch_assembly )
            .set { ch_for_polish }

        MINIMAP2_POLISH (
            ch_for_polish.map { meta, sr, lr, fasta -> tuple(meta, lr)  },
            ch_for_polish.map { meta, sr, lr, fasta -> fasta  },
            true,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_POLISH.out.versions.ifEmpty(null))

        SAMTOOLS_INDEX (
            MINIMAP2_POLISH.out.bam.dump(tag: 'samtools_sort')
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.ifEmpty(null))

        ch_for_polish    // tuple val(meta), val(reads), file(longreads), file(assembly)
            .join( MINIMAP2_POLISH.out.bam )    // tuple val(meta), file(bam)
            .join( SAMTOOLS_INDEX.out.bai )     // tuple  val(meta), file(bai)
            .join( ch_fast5 )             // tuple val(meta), file(fast5)
            .set { ch_for_nanopolish }          // tuple val(meta), val(reads), file(longreads), file(assembly), file(bam), file(bai), file(fast5)

        // TODO: 'nanopolish index' couldn't be tested. No fast5 provided in test datasets.
        NANOPOLISH (
            ch_for_nanopolish.dump(tag: 'into_nanopolish')
        )
        ch_versions = ch_versions.mix(NANOPOLISH.out.versions.ifEmpty(null))
    }

    //
    // MODULE: Medaka, polishes assembly - should take either miniasm, canu, or unicycler consensus sequence
    //
    if ( !params.skip_polish && params.assembly_type == 'long' && params.polish_method == 'medaka' ) {
        ch_for_assembly
            .join( ch_assembly )
            .map { meta, sr, lr, assembly -> tuple(meta, lr, assembly) }
            .set { ch_for_medaka }

        MEDAKA ( ch_for_medaka.dump(tag: 'into_medaka') )
        ch_versions = ch_versions.mix(MEDAKA.out.versions.ifEmpty(null))
    }
*/
    //
    // MODULE: Kraken2, QC for sample purity
    //
    ch_kraken_short_multiqc = Channel.empty()
    ch_kraken_long_multiqc  = Channel.empty()
    if ( !params.skip_kraken2 ) {
        KRAKEN2_DB_PREPARATION (
            kraken2db
        )
        KRAKEN2 (
            ch_for_kraken2_short.dump(tag: 'kraken2_short'),
            KRAKEN2_DB_PREPARATION.out.db.map { info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_short_multiqc = KRAKEN2.out.report
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.ifEmpty(null))

        KRAKEN2_LONG (
            ch_for_kraken2_long
                .map { meta, reads ->
                    info = [:]
                    info.id = meta.id
                    info.single_end = true
                    [ info, reads ]
                }
                .dump(tag: 'kraken2_long'),
            KRAKEN2_DB_PREPARATION.out.db.map { info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_long_multiqc = KRAKEN2_LONG.out.report
        ch_versions = ch_versions.mix(KRAKEN2_LONG.out.versions.ifEmpty(null))
    }

    //
    // MODULE: Kmerfinder, QC for sample purity
    //

    // TODO: add check contamination module // CALLIT PARSE_KMERFINDER
    // TODO: if not provided, download reference from kmerfinder results --> module FIND_DOWNLOAD_COMMON_REFFERENCE
    // TODO: Create kmerfinder mode for short and longreads
    // TODO: When no kmerfinder database is found, allow nf-core/bacass to download it
    if ( !params.skip_kmerfinder && params.kmerfinderdb ) {
        KMERFINDER (
            ch_for_assembly.map{ meta, sr, lr -> tuple( meta, sr) },    // [meta, reads]
            params.kmerfinderdb // path(kmerfinder database)
        )
        ch_kmerfinder_report  = KMERFINDER.out.report
        ch_versions           = ch_versions.mix( KMERFINDER.out.versions.ifEmpty(null) )

        KMERFINDER_SUMMARY (
            ch_kmerfinder_report.map{meta, report -> report }.collect()
        )
        ch_versions = ch_versions.mix( KMERFINDER_SUMMARY.out.versions.ifEmpty(null) )
    }
/*
    //
    // MODULE: QUAST, assembly QC
    //
    ch_assembly
        .collect{ it[1] }
        .map { consensus_collect -> tuple([id: "report"], consensus_collect) }
        .set { ch_to_quast }

    QUAST (
        ch_to_quast,
        [[:],[]],
        [[:],[]]
    )
    ch_quast_multiqc = QUAST.out.tsv
    ch_versions      = ch_versions.mix(QUAST.out.versions.ifEmpty(null))

    // Check assemblies that require further processing for gene annotation
    ch_assembly
        .branch{ meta, fasta ->
            gzip: fasta.name.endsWith('.gz')
            skip: true
        }
        .set{ ch_assembly_for_gunzip }

    //
    // MODULE: PROKKA, gene annotation
    //
    ch_prokka_txt_multiqc = Channel.empty()
    if ( !params.skip_annotation && params.annotation_tool == 'prokka' ) {
        // Uncompress assembly for annotation if necessary
        GUNZIP ( ch_assembly_for_gunzip.gzip )
        ch_to_prokka    = ch_assembly_for_gunzip.skip.mix( GUNZIP.out.gunzip )
        ch_versions     = ch_versions.mix( GUNZIP.out.versions.ifEmpty(null) )

        PROKKA (
            ch_to_prokka,
            [],
            []
        )
        ch_prokka_txt_multiqc   = PROKKA.out.txt.collect()
        ch_versions             = ch_versions.mix(PROKKA.out.versions.ifEmpty(null))
    }

    //
    // MODULE: BAKTA, gene annotation
    //
    ch_bakta_txt_multiqc = Channel.empty()
    if ( !params.skip_annotation && params.annotation_tool == 'bakta' ) {
        // Uncompress assembly for annotation if necessary
        GUNZIP ( ch_assembly_for_gunzip.gzip )
        ch_to_bakta     = ch_assembly_for_gunzip.skip.mix( GUNZIP.out.gunzip )
        ch_versions     = ch_versions.mix( GUNZIP.out.versions.ifEmpty(null) )

        BAKTA_DBDOWNLOAD_RUN (
            ch_to_bakta,
            params.baktadb,
            params.baktadb_download
        )
        ch_bakta_txt_multiqc    = BAKTA_DBDOWNLOAD_RUN.out.bakta_txt_multiqc.collect()
        ch_versions             = ch_versions.mix(BAKTA_DBDOWNLOAD_RUN.out.versions)
    }
    //
    // MODULE: DFAST, gene annotation
    //
    // TODO: "dfast_file_downloader.py --protein dfast --dbroot ." could be used in a separate process and the db could be forwarded
    if ( !params.skip_annotation && params.annotation_tool == 'dfast' ) {
        DFAST (
            ch_assembly,
            Channel.value(params.dfast_config ? file(params.dfast_config) : "")
        )
        ch_versions = ch_versions.mix(DFAST.out.versions.ifEmpty(null))
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc){
        workflow_summary        = WorkflowBacass.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary     = Channel.value(workflow_summary)
        methods_description     = WorkflowBacass.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description  = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_trim_json_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_kraken_short_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_kraken_long_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_quast_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_prokka_txt_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_bakta_txt_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_nanoplot_txt_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_porechop_log_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_pycoqc_multiqc.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
*/
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
