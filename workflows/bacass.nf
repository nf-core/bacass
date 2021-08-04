/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBacass.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { SKEWER                } from '../modules/local/skewer'                addParams( options: modules['skewer']            )
include { PORECHOP              } from '../modules/local/porechop'              addParams( options: modules['porechop']          )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC    } from '../modules/nf-core/modules/fastqc/main'    addParams( options: modules['fastqc']    )
include { NANOPLOT  } from '../modules/nf-core/modules/nanoplot/main'  addParams( options: modules['nanoplot']  )
include { PYCOQC    } from '../modules/nf-core/modules/pycoqc/main'    addParams( options: modules['pycoqc']    )
include { UNICYCLER } from '../modules/nf-core/modules/unicycler/main' addParams( options: modules['unicycler'] )
include { SPADES    } from '../modules/nf-core/modules/spades/main'    addParams( options: modules['spades']    )
include { MULTIQC   } from '../modules/nf-core/modules/multiqc/main'   addParams( options: multiqc_options      )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BACASS {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.shortreads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Skewer, trim and combine short read read-pairs per sample. Similar to nf-core vipr
    //
    SKEWER (
        INPUT_CHECK.out.shortreads
    )
    //ch_software_versions = ch_software_versions.mix(SKEWER.out.version.first().ifEmpty(null)) //TODO

    //
    // MODULE: Nanoplot, quality check for nanopore reads and Quality/Length Plots
    //
    NANOPLOT (
        INPUT_CHECK.out.longreads
    )
    ch_software_versions = ch_software_versions.mix(NANOPLOT.out.version.first().ifEmpty(null))

    //
    // MODULE: PYCOQC, quality check for nanopore reads and Quality/Length Plots
    //
    PYCOQC (
        INPUT_CHECK.out.fast5
    )
    ch_software_versions = ch_software_versions.mix(PYCOQC.out.version.first().ifEmpty(null))

    //
    // MODULE: PYCOQC, quality check for nanopore reads and Quality/Length Plots
    //

    // TODO: if ( params.assembly_type == 'hybrid' || params.assembly_type == 'long' && !('short' in params.assembly_type) )
    PORECHOP (
        INPUT_CHECK.out.longreads
    )
    ch_software_versions = ch_software_versions.mix(PORECHOP.out.version.first().ifEmpty(null))

    //
    // MODULE: Unicycler, nf-core module allows only short assembly assembly
    //
    // TODO: if ( params.assembler == 'unicycler' ) {

    UNICYCLER (
        SKEWER.out.reads
    )
    ch_software_versions = ch_software_versions.mix(UNICYCLER.out.version.first().ifEmpty(null))
/*
    //
    // MODULE: SPAdes, nf-core module allows only short read assembly
    //
    // TODO: if ( params.assembler == 'spades' ) {
    // TODO: error: java.lang.NullPointerException: Cannot get property 'id' on null object
    // TODO: probabaly: params.spades_hmm = false?
    params.spades_hmm = false
    SPADES (
        SKEWER.out.reads, params.spades_hmm
    )
    ch_software_versions = ch_software_versions.mix(SPADES.out.version.first().ifEmpty(null))
*/

    //
    // Join channels for unicycler, as trimming the files happens in two separate processes for paralellization of individual steps. As samples have the same sampleID, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
    // TODO:
/*
    if(params.assembly_type == 'hybrid'){
        SKEWER.out.reads
            .join(PORECHOP.out.reads)
            .dump(tag: 'unicycler')
            .set {ch_unicycler}
    } else if(params.assembly_type == 'short'){
        SKEWER.out.reads
            .map{ id,R1,R2 -> tuple(id,R1,R2,'NA') }
            .dump(tag: 'unicycler')
            .set {ch_unicycler}
    } else if(params.assembly_type == 'long'){
        ch_long_trimmed_unicycler
            .map{ id,lr -> tuple(id,'NA','NA',lr) }
            .dump(tag: 'unicycler')
            .set {ch_unicycler} //old channel name: ch_short_long_joint_unicycler
    }
*/

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBacass.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
