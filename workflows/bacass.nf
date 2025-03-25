/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Local to the pipeline
//
include { PYCOQC                    } from '../modules/local/pycoqc'
include { NANOPOLISH                } from '../modules/local/nanopolish'
include { MEDAKA                    } from '../modules/local/medaka'
include { KRAKEN2_DB_PREPARATION    } from '../modules/local/kraken2/db_preparation'
include { DFAST                     } from '../modules/local/dfast'
include { CUSTOM_MULTIQC            } from '../modules/local/custom/multiqc'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                } from '../modules/nf-core/fastqc'
include { CAT_FASTQ                             } from '../modules/nf-core/cat/fastq'
include { NANOPLOT                              } from '../modules/nf-core/nanoplot'
include { PORECHOP_PORECHOP                     } from '../modules/nf-core/porechop/porechop'
include { UNICYCLER                             } from '../modules/nf-core/unicycler'
include { CANU                                  } from '../modules/nf-core/canu'
include { MINIMAP2_ALIGN                        } from '../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_CONSENSUS  } from '../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_POLISH     } from '../modules/nf-core/minimap2/align'
include { MINIASM                               } from '../modules/nf-core/miniasm'
include { DRAGONFLYE                            } from '../modules/nf-core/dragonflye'
include { RACON                                 } from '../modules/nf-core/racon'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index'
include { KRAKEN2_KRAKEN2 as KRAKEN2            } from '../modules/nf-core/kraken2/kraken2'
include { KRAKEN2_KRAKEN2 as KRAKEN2_LONG       } from '../modules/nf-core/kraken2/kraken2'
include { QUAST                                 } from '../modules/nf-core/quast'
include { QUAST as QUAST_BYREFSEQID             } from '../modules/nf-core/quast'
include { GUNZIP                                } from '../modules/nf-core/gunzip'
include { PROKKA                                } from '../modules/nf-core/prokka'
include { FILTLONG                              } from '../modules/nf-core/filtlong'


//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_TRIM_FASTP_FASTQC               } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc'
include { KMERFINDER_SUMMARY_DOWNLOAD           } from '../subworkflows/local/kmerfinder_summary_download'
include { BAKTA_DBDOWNLOAD_RUN                  } from '../subworkflows/local/bakta_dbdownload_run'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_bacass_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db, params.dfast_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACASS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    def criteria = multiMapCriteria {
        meta, fastqs, long_fastq, fast5 ->
            shortreads: meta.single_end != 'NA' ? tuple(meta, fastqs) : null
            longreads: long_fastq       != 'NA' ? tuple(meta,long_fastq) : null
            fast5: fast5                != 'NA' ? tuple(meta, fast5) : null
    }
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    ch_samplesheet
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
    // MODULE: Concatenate FastQ files from same sample if required (shortreads)
    //
    ch_shortreads
        .branch{
            meta, fastqs ->
                single: fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_shortreads_fastqs }

        CAT_FASTQ (
            ch_shortreads_fastqs.multiple
        )
        .reads
        .mix( ch_shortreads_fastqs.single )
        .set { ch_shortreads_concat }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // SUBWORKFLOW: Short reads QC and trim adapters
    //
    ch_fastqc_raw_multiqc = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_fastp_json_multiqc = Channel.empty()
    if (params.assembly_type != 'long'){
        FASTQ_TRIM_FASTP_FASTQC (
        ch_shortreads_concat,
        [],
        params.save_trimmed_fail,
        [],
        params.skip_fastp,
        params.skip_fastqc
        )
        ch_fastqc_raw_multiqc   = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc  = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip
        ch_fastp_json_multiqc   = FASTQ_TRIM_FASTP_FASTQC.out.trim_json
        ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)
    }

    //
    // MODULE: Nanoplot, quality check for nanopore reads and Quality/Length Plots
    //
    NANOPLOT (
        ch_longreads
    )
    ch_nanoplot_txt_multiqc = NANOPLOT.out.txt
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)


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
        ch_versions       = ch_versions.mix(PYCOQC.out.versions)
    }

    //
    // MODULE: PORECHOP, quality check for nanopore reads and Quality/Length Plots
    //
    ch_porechop_log_multiqc = Channel.empty()
    if ((params.assembly_type == 'hybrid' || params.assembly_type == 'long' && !('short' in params.assembly_type)) && params.long_reads_filtering == 'porechop' ) {
        PORECHOP_PORECHOP (
            ch_longreads.dump(tag: 'longreads')
        )
        filtered_long_reads = PORECHOP_PORECHOP.out.reads
        ch_porechop_log_multiqc = PORECHOP_PORECHOP.out.log
        ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions )
    }

    //
    // MODULE: FILTLONG, filtering long reads by quality. It can take a set of long reads and produce a smaller, better subset.
    //
    ch_filtlong_multiqc = Channel.empty()
    //if ((params.assembly_type == 'hybrid' || params.assembly_type == 'long' && !('short' in params.assembly_type)) && params.long_reads_filtering == 'filtlong' ) {
    if ( !('short' in params.assembly_type) && params.long_reads_filtering == 'filtlong' ) { 
        if (params.assembly_type == 'hybrid') {
            ch_shortreads_4_filtlong = FASTQ_TRIM_FASTP_FASTQC.out.reads
        } else if ( params.assembly_type == 'long' ) {
            ch_shortreads_4_filtlong = ch_longreads.map{ meta, reads -> tuple(meta, []) }
        }

        FILTLONG (
            ch_shortreads_4_filtlong.join(ch_longreads)
        )

        filtered_long_reads = FILTLONG.out.reads
        ch_filtlong_multiqc = FILTLONG.out.log
        ch_versions       = ch_versions.mix(FILTLONG.out.versions)
    }


    //
    // Join channels for assemblers. As samples have the same meta data, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
    // Prepare channel for Kraken2
    //
    if(params.assembly_type == 'hybrid'){
        ch_for_kraken2_short    = FASTQ_TRIM_FASTP_FASTQC.out.reads
        ch_for_kraken2_long     = filtered_long_reads
        FASTQ_TRIM_FASTP_FASTQC.out.reads
            .dump(tag: 'fastp')
            .join(filtered_long_reads)
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
        ch_for_kraken2_long     = filtered_long_reads
        filtered_long_reads
            .dump(tag: 'filtered_long_reads')
            .map{ meta,lr -> tuple(meta,[],lr) }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
    }

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
        ch_versions = ch_versions.mix( UNICYCLER.out.versions )
    }


    //
    // MODULE: Canu, genome assembly, long reads
    //
    if ( params.assembler == 'canu' ) {
        CANU (
            ch_for_assembly.map { meta, reads, lr -> tuple( meta, lr ) },
            params.canu_mode,
            ch_for_assembly.map { meta, reads, lr -> meta.gsize }
        )
        ch_assembly = ch_assembly.mix( CANU.out.assembly.dump(tag: 'canu') )
        ch_versions = ch_versions.mix(CANU.out.versions)
    }

    //
    // MODULE: Miniasm, genome assembly, long reads
    //
    if ( params.assembly_type != 'short' && params.assembler == 'miniasm' ) {
        MINIMAP2_ALIGN (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            [[:],[]],
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        ch_for_assembly
            .join(MINIMAP2_ALIGN.out.paf)
            .map { meta, sr, lr, paf-> tuple(meta, lr, paf) }
            .set { ch_for_miniasm }

        MINIASM (
            ch_for_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions)

        MINIMAP2_CONSENSUS (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            MINIASM.out.assembly,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS.out.versions)

        ch_for_assembly
            .join(MINIASM.out.assembly)
            .join(MINIMAP2_CONSENSUS.out.paf)
            .map { meta, sr, lr, assembly, paf -> tuple(meta, lr, assembly, paf) }
            .set{ ch_for_racon }

        RACON (
            ch_for_racon
        )
        ch_assembly = ch_assembly.mix( RACON.out.improved_assembly.dump(tag: 'miniasm') )
        ch_versions = ch_versions.mix( RACON.out.versions )
    } else if (params.assembly_type == 'short' && params.assembler == 'miniasm') {
        exit("Selected assembler ${params.assembler} cannot run on short reads")
    }

    //
    // MODULE: Dragonflye, genome assembly of long reads. Moreover, it provides the option for polishing the draft genome using short reads when both short and long reads are available.
    //
    if( params.assembler == 'dragonflye' ){
        DRAGONFLYE(
            ch_for_assembly
        )
        ch_assembly = ch_assembly.mix( DRAGONFLYE.out.contigs.dump(tag: 'dragonflye') )
        ch_versions = ch_versions.mix( DRAGONFLYE.out.versions )
    }

    //
    // SUBWORKFLOW: Long reads polishing. Uses medaka or Nanopolish (this last requires Fast5 files available in input samplesheet).
    //
    if ( (params.assembly_type == 'long' && !params.skip_polish) || ( params.assembly_type != 'short' && params.polish_method) ){
        // Set channel for polishing long reads
        ch_for_assembly
            .join( ch_assembly )
            .map { meta, sr, lr, fasta -> tuple(meta, lr, fasta) }
            .set { ch_polish_long } // channel: [ val(meta), path(lr), path(fasta) ]
        if (params.polish_method == 'medaka'){
            //
            // MODULE: Medaka, polishes assembly - should take either miniasm, canu, or unicycler consensus sequence
            //
            MEDAKA ( ch_polish_long )
            ch_assembly = MEDAKA.out.assembly
            ch_versions = ch_versions.mix(MEDAKA.out.versions)
        } else if (params.polish_method == 'nanopolish') {
            //
            // MODULE: Nanopolish, polishes assembly using FAST5 files
            //
            if (!ch_fast5){
                log.error "ERROR: FAST5 files are required for Nanopolish but none were provided. Please supply FAST5 files or choose another polishing method. Available options are: medaka, nanopolish"
            } else {
                //
                // MODULE: Minimap2 polish
                //
                MINIMAP2_POLISH (
                    ch_polish_long.map { meta, lr, fasta -> tuple(meta, lr) },
                    ch_polish_long.map { meta, lr, fasta -> tuple(meta, fasta) },
                    true,
                    false,
                    false
                )
                ch_versions = ch_versions.mix(MINIMAP2_POLISH.out.versions)
                //
                // MODULE: Samtools index
                //
                SAMTOOLS_INDEX (
                    MINIMAP2_POLISH.out.bam.dump(tag: 'samtools_sort')
                )
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
                //
                // MODULE: Nanopolish
                //
                ch_polish_long                         // tuple val(meta), val(reads), file(longreads), file(assembly)
                    .join( MINIMAP2_POLISH.out.bam )  // tuple val(meta), file(bam)
                    .join( SAMTOOLS_INDEX.out.bai )   // tuple  val(meta), file(bai)
                    .join( ch_fast5 )                 // tuple val(meta), file(fast5)
                    .set { ch_for_nanopolish }        // tuple val(meta), val(reads), file(longreads), file(assembly), file(bam), file(bai), file(fast5)
                // TODO: 'nanopolish index' couldn't be tested. No fast5 provided in test datasets.
                NANOPOLISH (
                    ch_for_nanopolish.dump(tag: 'into_nanopolish')
                )
                ch_assembly = NANOPOLISH.out.assembly
                ch_versions = ch_versions.mix( NANOPOLISH.out.versions )
            }
        }
    }

    //
    // MODULE: Kraken2, QC for sample purity
    //
    ch_kraken_short_multiqc = Channel.empty()
    ch_kraken_long_multiqc  = Channel.empty()
    if ( !params.skip_kraken2 ) {
        KRAKEN2_DB_PREPARATION (
            params.kraken2db
        )
        ch_versions = ch_versions.mix(KRAKEN2_DB_PREPARATION.out.versions)
        KRAKEN2 (
            ch_for_kraken2_short.dump(tag: 'kraken2_short'),
            KRAKEN2_DB_PREPARATION.out.db.map { info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_short_multiqc = KRAKEN2.out.report
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

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
        ch_versions = ch_versions.mix(KRAKEN2_LONG.out.versions)
    }

    //
    // SUBWORKFLOW: Kmerfinder, QC for sample purity.
    //
    // Executes both kmerfinder and classifies samples by their reference genome (all this through the KMERFINDER_SUMMARY_DOWNLOAD()).

    ch_kmerfinder_multiqc = Channel.empty()
    if (!params.skip_kmerfinder) {
        // Set kmerfinder channel based on assembly type
        if( params.assembly_type == 'short' || params.assembly_type == 'hybrid' ) {
            ch_for_kmerfinder = FASTQ_TRIM_FASTP_FASTQC.out.reads
        } else if ( params.assembly_type == 'long' ) {
            ch_for_kmerfinder = filtered_long_reads
        }
        // RUN kmerfinder subworkflow
        KMERFINDER_SUMMARY_DOWNLOAD (
            ch_for_kmerfinder,
            ch_assembly
        )
        ch_kmerfinder_multiqc   = KMERFINDER_SUMMARY_DOWNLOAD.out.summary_yaml
        ch_consensus_byrefseq   = KMERFINDER_SUMMARY_DOWNLOAD.out.consensus_byrefseq
        ch_versions             = ch_versions.mix(KMERFINDER_SUMMARY_DOWNLOAD.out.versions)

        // Set channel to perform by refseq QUAST based on reference genome identified with KMERFINDER.
        ch_consensus_byrefseq
            .map {
                refmeta, meta, consensus, ref_fna, ref_gff ->
                    return tuple(refmeta, consensus.flatten(), ref_fna, ref_gff)
            }
            .set { ch_to_quast_byrefseq }
    }

    //
    // MODULE: QUAST, assembly QC
    //
    ch_assembly
        .collect{it[1]}
        .map{ consensus -> tuple([id:'report'], consensus) }
        .set{ ch_to_quast }

    if(params.skip_kmerfinder){
        QUAST(
            ch_to_quast,
            params.reference_fasta ?: [[:],[]],
            params.reference_gff ?: [[:],[]]
        )
        ch_quast_multiqc = QUAST.out.results
    } else if (!params.skip_kmerfinder) {
        // Quast runs twice if kmerfinder is allowed.
        // This approach allow Quast to calculate relevant parameters such as genome fraction based on a reference genome.
        QUAST(
            ch_to_quast,
            [[:],[]],
            [[:],[]]
        )
        QUAST_BYREFSEQID(
            ch_to_quast_byrefseq.map{ refmeta, consensus, ref_fasta, ref_gff -> tuple( refmeta, consensus)},
            ch_to_quast_byrefseq.map{ refmeta, consensus, ref_fasta, ref_gff -> tuple( refmeta, ref_fasta)},
            ch_to_quast_byrefseq.map{ refmeta, consensus, ref_fasta, ref_gff -> tuple( refmeta, ref_gff)}
        )
        ch_quast_multiqc = QUAST_BYREFSEQID.out.results
        ch_versions      = ch_versions.mix(QUAST_BYREFSEQID.out.versions)
    }
    ch_versions = ch_versions.mix(QUAST.out.versions)

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
        ch_versions     = ch_versions.mix( GUNZIP.out.versions )

        PROKKA (
            ch_to_prokka.filter{ meta, fasta -> !fasta.isEmpty() },
            [],
            []
        )
        ch_prokka_txt_multiqc   = PROKKA.out.txt.map{ meta, prokka_txt -> [ prokka_txt ]}
        ch_versions             = ch_versions.mix(PROKKA.out.versions)
    }

    //
    // MODULE: BAKTA, gene annotation
    //
    ch_bakta_txt_multiqc = Channel.empty()
    if ( !params.skip_annotation && params.annotation_tool == 'bakta' ) {
        // Uncompress assembly for annotation if necessary
        GUNZIP ( ch_assembly_for_gunzip.gzip )
        ch_to_bakta     = ch_assembly_for_gunzip.skip.mix( GUNZIP.out.gunzip )
        ch_versions     = ch_versions.mix( GUNZIP.out.versions )

        BAKTA_DBDOWNLOAD_RUN (
            ch_to_bakta.filter{ meta, fasta -> !fasta.isEmpty() },
            params.baktadb,
            params.baktadb_download
        )
        ch_bakta_txt_multiqc    = BAKTA_DBDOWNLOAD_RUN.out.bakta_txt_multiqc.map{ meta, bakta_txt -> [ bakta_txt ]}
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
        ch_versions = ch_versions.mix(DFAST.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'bacass_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = !params.skip_kmerfinder && params.assembly_type ? Channel.fromPath("$projectDir/assets/multiqc_config_${params.assembly_type}.yml", checkIfExists: true) : Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? Channel.fromPath(params.multiqc_methods_description, checkIfExists: true) : Channel.fromPath("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    CUSTOM_MULTIQC (
        ch_multiqc_config.ifEmpty([]),
        ch_multiqc_custom_config.ifEmpty([]),
        ch_multiqc_logo.ifEmpty([]),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        ch_multiqc_custom_methods_description.ifEmpty([]),
        ch_collated_versions.ifEmpty([]),
        ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
        ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
        ch_fastp_json_multiqc.collect{it[1]}.ifEmpty([]),
        ch_nanoplot_txt_multiqc.collect{it[1]}.ifEmpty([]),
        ch_porechop_log_multiqc.collect{it[1]}.ifEmpty([]),
        ch_pycoqc_multiqc.collect{it[1]}.ifEmpty([]),
        ch_kraken_short_multiqc.collect{it[1]}.ifEmpty([]),
        ch_kraken_long_multiqc.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect{it[1]}.ifEmpty([]),
        ch_prokka_txt_multiqc.collect().ifEmpty([]),
        ch_bakta_txt_multiqc.collect().ifEmpty([]),
        ch_kmerfinder_multiqc.collectFile(name: 'multiqc_kmerfinder.yaml').ifEmpty([]),
    )
    multiqc_report = CUSTOM_MULTIQC.out.report.toList()

    emit:
    multiqc_report = CUSTOM_MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                        // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
