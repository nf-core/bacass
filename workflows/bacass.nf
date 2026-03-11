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
include { CAT_FASTQ as CAT_FASTQ_SHORT          } from '../modules/nf-core/cat/fastq'
include { CAT_FASTQ as CAT_FASTQ_LONG           } from '../modules/nf-core/cat/fastq'
include { PORECHOP_PORECHOP                     } from '../modules/nf-core/porechop/porechop'
include { AUTOCYCLER_SUBSAMPLE                  } from '../modules/nf-core/autocycler/subsample/main'
include { UNICYCLER                             } from '../modules/nf-core/unicycler'
include { MEGAHIT                               } from '../modules/nf-core/megahit/main'
include { CANU                                  } from '../modules/nf-core/canu'
include { MINIMAP2_ALIGN                        } from '../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_CONSENSUS  } from '../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_POLISH     } from '../modules/nf-core/minimap2/align'
include { MINIASM                               } from '../modules/nf-core/miniasm'
include { DRAGONFLYE                            } from '../modules/nf-core/dragonflye'
include { RAVEN                                 } from '../modules/nf-core/raven/main'
include { FLYE                                  } from '../modules/nf-core/flye/main'
include { RACON                                 } from '../modules/nf-core/racon'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index'
include { KRAKEN2_KRAKEN2 as KRAKEN2            } from '../modules/nf-core/kraken2/kraken2'
include { KRAKEN2_KRAKEN2 as KRAKEN2_LONG       } from '../modules/nf-core/kraken2/kraken2'
include { QUAST                                 } from '../modules/nf-core/quast'
include { QUAST as QUAST_BYREFSEQID             } from '../modules/nf-core/quast'
include { QUAST as QUAST_BYSAMPLE               } from '../modules/nf-core/quast'
include { BUSCO_BUSCO                           } from '../modules/nf-core/busco/busco/main'
include { GUNZIP                                } from '../modules/nf-core/gunzip'
include { PROKKA                                } from '../modules/nf-core/prokka'
include { FILTLONG                              } from '../modules/nf-core/filtlong'
include { LIFTOFF                               } from '../modules/nf-core/liftoff'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//

include { FASTQ_TRIM_FASTP_FASTQC               } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { QC_NANOPLOT_TOULLIGQC                 } from '../subworkflows/local/qc_nanoplot_toulliqc'
include { KMERFINDER_SUMMARY_DOWNLOAD           } from '../subworkflows/local/kmerfinder_summary_download'
include { FASTA_CONSENSUS_AUTOCYCLER            } from '../subworkflows/nf-core/fasta_consensus_autocycler/main'
include { BAKTA_DBDOWNLOAD_RUN                  } from '../subworkflows/local/bakta_dbdownload_run'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_bacass_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACASS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db, params.dfast_config, params.reference_fasta, params.reference_gff ]
    checkPathParamList.each { param -> if (param) { file(param, checkIfExists: true) } }

    if (params.reference_fasta) {
        reference_fasta = file(params.reference_fasta, type: 'file')
    }
    if (params.reference_gff) {
        reference_gff = file(params.reference_gff, type: 'file')
    }


    ch_versions = channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    def criteria = multiMapCriteria {
        meta, fastqs, long_fastq, fast5 ->
            shortreads: fastqs[0]       != 'NA' ? tuple(meta, fastqs)    : tuple(meta, [])
            longreads:  long_fastq      != 'NA' ? tuple(meta,long_fastq) : tuple(meta, [])
            fast5:      fast5           != 'NA' ? tuple(meta, fast5)     : tuple(meta, [])
    }
    ch_proteins = params.prokka_proteins ? channel.fromPath(params.prokka_proteins, checkIfExists: true)  : []
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    ch_samplesheet
        .map { meta, fastqs, long_fastq, fast5  ->
            def new_meta = meta + [id: meta.sample] // add "meta.id" !
            new_meta.subsample = false
            return [ new_meta, fastqs, long_fastq, fast5  ] }
        .multiMap (criteria)
        .set { ch_input }

    // reconfigure channels
    ch_input
        .shortreads
        .filter{ meta, data -> data }
        .set { ch_shortreads }
    ch_input
        .longreads
        .filter{ meta, data -> data }
        .set { ch_longreads }
    ch_input
        .fast5
        .filter{ meta, data -> data }
        .set { ch_fast5 }

    //
    // Short read preprocessing and QC
    //
    ch_short_preprocessed  = channel.empty()
    ch_fastqc_raw_multiqc  = channel.empty()
    ch_fastqc_trim_multiqc = channel.empty()
    ch_fastp_json_multiqc  = channel.empty()
    if (params.assembly_type in ['short', 'hybrid']) {
        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        ch_shortreads
            .branch{
                meta, fastqs ->
                    single: meta.single_end ? fastqs.size() == 1 : fastqs.size() == 2
                    multiple: meta.single_end ? fastqs.size() > 1 : fastqs.size() > 2
            }
            .set { ch_shortreads_fastqs }
        CAT_FASTQ_SHORT (
            ch_shortreads_fastqs.multiple
        )
        ch_shortreads_concat = CAT_FASTQ_SHORT.out.reads
            .mix( ch_shortreads_fastqs.single )

        //
        // SUBWORKFLOW: Short reads QC and trim adapters
        //
        FASTQ_TRIM_FASTP_FASTQC (
            ch_shortreads_concat.map{ meta, reads -> [meta, reads, []] }, //add an empty field for adapters
            params.save_trimmed_fail,
            [],
            params.discard_trimmed_pass,
            params.skip_fastp,
            params.skip_fastqc
        )
        ch_short_preprocessed   = FASTQ_TRIM_FASTP_FASTQC.out.reads
        ch_fastqc_raw_multiqc   = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc  = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip
        ch_fastp_json_multiqc   = FASTQ_TRIM_FASTP_FASTQC.out.trim_json
        ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)
    }

    //
    // Long read preprocessing and QC
    //
    ch_pycoqc_multiqc       = channel.empty()
    ch_nanoplot_txt_multiqc = channel.empty()
    ch_porechop_log_multiqc = channel.empty()
    ch_filtlong_log_multiqc = channel.empty()
    ch_longreads_filtered   = channel.empty()
    if (params.assembly_type in ['long', 'hybrid']) {
        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        ch_longreads
            .map {
                meta, long_fastq ->
                    // Force single_end=true for long reads
                    // Create a copy of meta to avoid interference with short reads meta (when hybrid mode is activated)
                    def new_meta = meta + [single_end: true]  // Force single_end
                    return [ new_meta, long_fastq ]
            }
            .branch{
                _meta, long_fastqs ->
                    single: long_fastqs.size() == 1
                    multiple: long_fastqs.size() > 1
            }
            .set { ch_longreads_fastqs }
        CAT_FASTQ_LONG (
            ch_longreads_fastqs.multiple
        )
        ch_longreads_concat = CAT_FASTQ_LONG.out.reads
            .mix( ch_longreads_fastqs.single )

        //
        // SUBWORKFLOW: quality check for nanopore reads with Nanoplot and ToulligQC
        //
        QC_NANOPLOT_TOULLIGQC (
            ch_longreads_concat,
            params.skip_nanoplot,  // skip the nanoplot qc
            params.skip_toulligqc  // skip the toulligqc
        )
        ch_nanoplot_txt_multiqc = QC_NANOPLOT_TOULLIGQC.out.nanoplot_txt
        ch_versions = ch_versions.mix(QC_NANOPLOT_TOULLIGQC.out.nanoplot_version)
        ch_versions = ch_versions.mix(QC_NANOPLOT_TOULLIGQC.out.toulligqc_version)

        //
        // MODULE: PYCOQC, quality check for nanopore reads and Quality/Length Plots
        //
        // TODO: Couldn't be tested. No configuration test available (lack of fast5 file or params.skip_pycoqc=false).
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
        if ( params.long_reads_filtering == 'porechop' ) {
            PORECHOP_PORECHOP (
                ch_longreads_concat.dump(tag: 'longreads')
            )
            ch_longreads_filtered   = PORECHOP_PORECHOP.out.reads
            ch_porechop_log_multiqc = PORECHOP_PORECHOP.out.log
            ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions )
        }

        //
        // MODULE: FILTLONG, filtering long reads by quality. It can take a set of long reads and produce a smaller, better subset.
        //
        if ( params.long_reads_filtering == 'filtlong' ) {
            ch_shortreads_for_filtlong = channel.empty()
            if (params.assembly_type == 'hybrid') {
                ch_shortreads_for_filtlong = ch_short_preprocessed.join(ch_longreads_concat)   //tuple val(meta), file(sr), file(lr)
            } else if ( params.assembly_type == 'long' ) {
                ch_shortreads_for_filtlong = ch_longreads_concat.map{ meta, lr -> tuple(meta, [], lr ) }
            }

            FILTLONG (
                ch_shortreads_for_filtlong
            )

            ch_longreads_filtered   = FILTLONG.out.reads
            ch_filtlong_log_multiqc = FILTLONG.out.log
            ch_versions       = ch_versions.mix(FILTLONG.out.versions)
        }
    }


    //
    // Join channels for assemblers. As samples have the same meta data, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
    // Prepare channel for Kraken2
    //
    if(params.assembly_type == 'hybrid'){
        ch_for_kraken2_short    = ch_short_preprocessed
        ch_for_kraken2_long     = ch_longreads_filtered
        ch_short_preprocessed
            .dump(tag: 'fastp')
            .cross(ch_longreads_filtered) { it -> it[0].sample }
            .map { short_tuple, long_tuple ->
                def meta_short = short_tuple[0]
                def short_reads = short_tuple[1]
                def long_reads = long_tuple[1]
                [meta_short, short_reads, long_reads]
            }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
        ch_for_assembly.ifEmpty{
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is nothing to assemble with these settings.\n" +
                "  Please verify that samples have short and long reads.\n"
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            error(error_string) }
    } else if ( params.assembly_type == 'short' ) {
        ch_for_kraken2_short    = ch_short_preprocessed
        ch_for_kraken2_long     = channel.empty()
        ch_short_preprocessed
            .dump(tag: 'fastp')
            .map{ meta,reads -> tuple(meta,reads,[]) }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
        ch_for_assembly.ifEmpty{
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is nothing to assemble with these settings.\n" +
                "  Please verify that samples have short reads.\n"
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            error(error_string) }
    } else if ( params.assembly_type == 'long' ) {
        ch_for_kraken2_short    = channel.empty()
        ch_for_kraken2_long     = ch_longreads_filtered
        ch_longreads_filtered
            .dump(tag: 'ch_longreads_filtered')
            .map{ meta,lr -> tuple(meta,[],lr) }
            .dump(tag: 'ch_for_assembly')
            .set { ch_for_assembly }
        ch_for_assembly.ifEmpty{
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is nothing to assemble with these settings.\n" +
                "  Please verify that samples have long reads.\n"
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            error(error_string) }
    }

    //
    // MODULE: Autocycler, subset long reads for multiple assemblies per sample.
    //
    if( params.assembler.tokenize(",").contains("autocycler") ) {
        // subsample and transpose to one subset per channel entry
        AUTOCYCLER_SUBSAMPLE (
            ch_for_assembly.map{ meta, _short_reads, long_reads -> [meta, long_reads] },
            ch_for_assembly.map { meta, _reads, _lr -> meta.gsize }
        )
        AUTOCYCLER_SUBSAMPLE.out.subsampled_reads // transpose to [ meta, fasta-subset1, fasta-subset2, ... ]
            .transpose() // transpose to [ meta, fasta ]
            .map{ meta, long_reads ->
                def new_meta = meta.clone()
                new_meta.subsample = long_reads.getBaseName() -'.fastq'
                [ new_meta, [], long_reads ]
            }
            .set{ ch_subsets }
        // data subsets are mixed into ch_for_assembly
        ch_for_assembly = ch_for_assembly.mix(ch_subsets)
    }

    //
    // ASSEMBLY: Unicycler, Megahit, Canu, Miniasm, Dragonflye, Raven, Flye, Autocycler
    //
    ch_assembly = channel.empty()

    //
    // MODULE: Unicycler, genome assembly, nf-core module allows only short, long and hybrid assembly
    //
    if ( params.assembler.tokenize(",").contains("unicycler") ) {
        ch_for_assembly
            .filter{ meta, sr, lr -> !meta.subsample } // subsamples are not entering. i.e. anything with "meta.subsample"
            .map{ meta, sr, lr ->
                def new_meta = meta.clone()
                new_meta.assembler = "unicycler"
                new_meta.id = meta.id + "-unicycler"
                [ new_meta, sr, lr ]
            }
            .set { ch_for_assembly_uniycler }
        UNICYCLER (
            ch_for_assembly_uniycler
        )
        ch_assembly = ch_assembly.mix( UNICYCLER.out.scaffolds.dump(tag: 'unicycler') )
        ch_versions = ch_versions.mix( UNICYCLER.out.versions )
    }

    //
    // MODULE: MEGAHIT, genome assembly of short reads
    //
    if ( params.assembler.tokenize(",").contains("megahit") ) {
        ch_for_assembly
            .map { meta, short_reads, _long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "megahit"
                new_meta.id = meta.id + "-megahit"
                def reads1 = meta.single_end ? short_reads : short_reads[0]
                def reads2 = meta.single_end ? [] : short_reads[1]
                [ new_meta, reads1, reads2 ]
            }
            .filter { meta, reads1, _reads2 -> !meta.subsample && reads1 } // only non-subsampled samples with short reads
            .set { ch_for_assembly_megahit }

        MEGAHIT (
            ch_for_assembly_megahit
        )
        ch_assembly = ch_assembly.mix( MEGAHIT.out.contigs.dump(tag: 'megahit') )
        ch_versions = ch_versions.mix( MEGAHIT.out.versions_megahit )
    }


    //
    // MODULE: Canu, genome assembly, long reads
    //
    if ( params.assembler.tokenize(",").contains("canu") || ( params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("canu") ) ) {
        ch_for_assembly
            .map{ meta, _short_reads, long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "canu"
                new_meta.id = meta.subsample ? "${meta.id}-${meta.subsample}-canu" : meta.id + "-canu"
                [ new_meta, long_reads ]
            }
            .filter { meta, lr ->
                params.assembler.tokenize(",").contains("autocycler") && !params.assembler.tokenize(",").contains("canu") ? meta.subsample : // if with autocycler and not canu accept only subsamples
                    !params.assembler.tokenize(",").contains("autocycler") && params.assembler.tokenize(",").contains("canu") ? !meta.subsample : true // if without autocycler and with canu reject subsample
            }
            .set { ch_for_assembly_canu }
        CANU (
            ch_for_assembly_canu,
            params.canu_mode,
            ch_for_assembly_canu.map { meta, _lr -> meta.gsize }
        )
        ch_assembly = ch_assembly.mix( CANU.out.assembly.dump(tag: 'canu') )
        ch_versions = ch_versions.mix(CANU.out.versions)
    }

    //
    // MODULE: Miniasm, genome assembly, long reads
    //
    if ( params.assembler.tokenize(",").contains("miniasm") || ( params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("miniasm") ) ) {
        ch_for_assembly
            .map{ meta, _short_reads, long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "miniasm"
                new_meta.id = meta.subsample ? "${meta.id}-${meta.subsample}-miniasm" : meta.id + "-miniasm"
                [ new_meta, long_reads ]
            }
            .filter { meta, lr ->
                params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("miniasm") && params.assembler.tokenize(",").contains("miniasm") ? true : // if with autocycler and miniasm accept all data sets
                    params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("miniasm") && !params.assembler.tokenize(",").contains("miniasm") ? meta.subsample : // if with autocycler and not miniasm accept only subsamples
                    ( !params.assembler.tokenize(",").contains("autocycler") || !params.autocycler_assemblers.tokenize(",").contains("miniasm") ) && params.assembler.tokenize(",").contains("miniasm") ? !meta.subsample : false // if without autocycler and with miniasm reject subsample
            }
            .set { ch_for_assembly_miniasm }

        MINIMAP2_ALIGN (
            ch_for_assembly_miniasm,
            [[:],[]],
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        ch_for_assembly_miniasm
            .join(MINIMAP2_ALIGN.out.paf)
            .set { ch_for_miniasm }

        MINIASM (
            ch_for_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions)

        MINIMAP2_CONSENSUS (
            ch_for_assembly_miniasm,
            MINIASM.out.assembly,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS.out.versions)

        ch_for_assembly_miniasm
            .join(MINIASM.out.assembly)
            .join(MINIMAP2_CONSENSUS.out.paf)
            .set{ ch_for_racon }

        RACON (
            ch_for_racon
        )
        ch_assembly = ch_assembly.mix( RACON.out.improved_assembly.dump(tag: 'miniasm') )
        ch_versions = ch_versions.mix( RACON.out.versions )
    }

    //
    // MODULE: Dragonflye, genome assembly of long reads. Moreover, it provides the option for polishing the draft genome using short reads when both short and long reads are available.
    //
    if( params.assembler.tokenize(",").contains("dragonflye") ){
        ch_for_assembly
            .map{ meta, short_reads, long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "dragonflye"
                new_meta.id = meta.id + "-dragonflye"
                [ new_meta, short_reads, long_reads ]
            }
            .filter{ meta, sr, lr -> !meta.subsample } // subsamples are not entering. i.e. anything with "meta.subsample"
            .set { ch_for_assembly_dragonflye }

        DRAGONFLYE(
            ch_for_assembly_dragonflye
        )
        ch_assembly = ch_assembly.mix( DRAGONFLYE.out.contigs.dump(tag: 'dragonflye') )
        ch_versions = ch_versions.mix( DRAGONFLYE.out.versions )
    }

    //
    // MODULE: Raven, genome assembly of long reads.
    //
    if ( params.assembler.tokenize(",").contains("raven") || ( params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("raven") ) ) {
        ch_for_assembly
            .map{ meta, _short_reads, long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "raven"
                new_meta.id = meta.subsample ? "${meta.id}-${meta.subsample}-raven" : meta.id + "-raven"
                [ new_meta, long_reads ]
            }
            .filter { meta, lr ->
                params.assembler.tokenize(",").contains("autocycler") && !params.assembler.tokenize(",").contains("raven") ? meta.subsample : // if with autocycler and not raven accept only subsamples
                    !params.assembler.tokenize(",").contains("autocycler") && params.assembler.tokenize(",").contains("raven") ? !meta.subsample : true // if without autocycler and with raven reject subsample
            }
            .set { ch_for_assembly_raven }

        RAVEN (
            ch_for_assembly_raven
        )
        ch_assembly = ch_assembly.mix( RAVEN.out.fasta )
    }

    //
    // MODULE: Flye, genome assembly of long reads.
    //
    if ( params.assembler.tokenize(",").contains("flye") || ( params.assembler.tokenize(",").contains("autocycler") && params.autocycler_assemblers.tokenize(",").contains("flye") ) ) {
        ch_for_assembly
            .map{ meta, _short_reads, long_reads ->
                def new_meta = meta.clone()
                new_meta.assembler = "flye"
                new_meta.id = meta.subsample ? "${meta.id}-${meta.subsample}-flye" : meta.id + "-flye"
                [ new_meta, long_reads ]
            }
            .filter { meta, lr ->
                params.assembler.tokenize(",").contains("autocycler") && !params.assembler.tokenize(",").contains("flye") ? meta.subsample : // if with autocycler and not flye accept only subsamples
                    !params.assembler.tokenize(",").contains("autocycler") && params.assembler.tokenize(",").contains("flye") ? !meta.subsample : true // if without autocycler and with flye reject subsample
            }
            .set { ch_for_assembly_flye }
        FLYE (
            ch_for_assembly_flye,
            params.flye_mode
        )
        ch_assembly = ch_assembly.mix( FLYE.out.fasta )
    }

    //
    // SUBWORKFLOW: Autocycler, combine genome assembly of long reads.
    //
    if( params.assembler.tokenize(",").contains("autocycler") ){
        ch_assembly
            .filter{ meta, assembly -> meta.subsample } // only assemblies of subsamples pass, i.e. anything with "meta.subsample"
            .map{ meta, assembly ->
                def new_meta = meta.clone()
                new_meta.remove("subsample")
                new_meta.assembler = "autocycler"
                new_meta.id = meta.sample + "-autocycler"
                [ new_meta, assembly ]
            }
            .filter{ meta, assembly -> assembly.countLines() > 1 } // keep only non-empty assembly files
            .groupTuple() // group to "[ val(meta), [ fasta, fasta, ... ] ]"
            .set { ch_assembly_autocycler }

        FASTA_CONSENSUS_AUTOCYCLER (
                ch_assembly_autocycler
            )
        // combine assemblies with autocycler combined assembly
        ch_assembly
            .mix( FASTA_CONSENSUS_AUTOCYCLER.out.consensus_assembly )
            .set { ch_assembly }
    }

    // clean assemblies from subsamples
    ch_assembly
        .filter { meta, assembly -> !meta.subsample } // omit subsample assemblies
        .set { ch_assembly }

    //
    // SUBWORKFLOW: Long reads polishing. Uses medaka or Nanopolish (this last requires Fast5 files available in input samplesheet).
    //
    if ( (params.assembly_type == 'long' && !params.skip_polish) || ( params.assembly_type != 'short' && params.polish_method) ){
        // Set channel for polishing long reads
        ch_for_assembly
            .filter { meta, sr, lr -> !meta.subsample } // remove any subsamples
            .cross(ch_assembly) { it -> it[0].sample } // merge by meta.sample -> [[ meta, sr, lr ],[ meta, assembly ]]
            .map { for_assembly,assembly ->
                def meta_assembly  = assembly[0]
                def long_reads     = for_assembly[2]
                def fasta_assembly = assembly[1]
                [ meta_assembly, long_reads, fasta_assembly ] }
            .set { ch_polish_long } // channel: [ val(meta), path(lr), path(fasta) ]

        if (params.polish_method == 'medaka'){
            ch_polish_long
                .map{ meta, lr, assembly ->
                    def new_meta = meta.clone()
                    new_meta.polish = "medaka"
                    new_meta.id = meta.id + "-medaka"
                    [ new_meta, lr, assembly ]
                }
                .set { ch_polish_long_medaka }

            //
            // MODULE: Medaka, polishes assembly - should take either miniasm, canu, or unicycler consensus sequence
            //
            MEDAKA ( ch_polish_long_medaka )
            ch_assembly = MEDAKA.out.assembly
            ch_versions = ch_versions.mix(MEDAKA.out.versions)
        } else if (params.polish_method == 'nanopolish') {
            ch_polish_long
                .map{ meta, lr, assembly ->
                    def new_meta = meta.clone()
                    new_meta.polish = "nanopolish"
                    new_meta.id = meta.id + "-nanopolish"
                    [ new_meta, lr, assembly ]
                }
                .set { ch_polish_long_nanopolish }
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
                    ch_polish_long_nanopolish.map { meta, lr, _fasta -> tuple(meta, lr) },
                    ch_polish_long_nanopolish.map { meta, _lr, fasta -> tuple(meta, fasta) },
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
                ch_polish_long_nanopolish                     // tuple val(meta), file(longreads), file(assembly)
                    .join( MINIMAP2_POLISH.out.bam )          // tuple val(meta), file(bam)
                    .join( SAMTOOLS_INDEX.out.bai )           // tuple val(meta), file(bai)
                    .cross( ch_fast5 ) { it -> it[0].sample } // tuple val(meta), file(fast5) // meta differs here and needs a join on meta.sample!
                    .map { input_tuple, fast5_tuple ->
                        def meta       = input_tuple[0]
                        def long_reads = input_tuple[1]
                        def assembly   = input_tuple[2]
                        def bam        = input_tuple[3]
                        def bai        = input_tuple[4]
                        def fast5      = fast5_tuple[1]
                        [meta, long_reads, assembly, bam, bai, fast5]
                    }
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
    ch_kraken_short_multiqc = channel.empty()
    ch_kraken_long_multiqc  = channel.empty()
    if ( !params.skip_kraken2 ) {
        KRAKEN2_DB_PREPARATION (
            params.kraken2db
        )
        ch_versions = ch_versions.mix(KRAKEN2_DB_PREPARATION.out.versions)
        KRAKEN2 (
            ch_for_kraken2_short.dump(tag: 'kraken2_short'),
            KRAKEN2_DB_PREPARATION.out.db.map { _info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_short_multiqc = KRAKEN2.out.report

        KRAKEN2_LONG (
            ch_for_kraken2_long
                .map { meta, reads ->
                    def info = [:]
                    info.id = meta.id
                    info.single_end = true
                    [ info, reads ]
                }
                .dump(tag: 'kraken2_long'),
            KRAKEN2_DB_PREPARATION.out.db.map { _info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_long_multiqc = KRAKEN2_LONG.out.report
    }

    //
    // SUBWORKFLOW: Kmerfinder, QC for sample purity.
    //
    // Executes both kmerfinder and classifies samples by their reference genome (all this through the KMERFINDER_SUMMARY_DOWNLOAD()).

    ch_kmerfinder_multiqc = channel.empty()
    if (!params.skip_kmerfinder) {
        // Set kmerfinder channel based on assembly type
        if( params.assembly_type == 'short' || params.assembly_type == 'hybrid' ) {
            ch_for_kmerfinder = ch_short_preprocessed
        } else if ( params.assembly_type == 'long' ) {
            ch_for_kmerfinder = ch_longreads_filtered
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
                refmeta, _meta, consensus, ref_fna, ref_gff ->
                    return tuple(refmeta, consensus.flatten(), ref_fna, ref_gff)
            }
            .set { ch_to_quast_byrefseq }
    }

    //
    // MODULE: QUAST, assembly QC
    //
    ch_assembly
        .collect{it -> it[1]}
        .map{ consensus -> tuple([id:'report'], consensus) }
        .set{ ch_to_quast }

    // First, check if there are multiple distinct samples
    ch_assembly
        .map { meta, _consensus -> meta.sample } // Extract sample value
        .unique()                                // Get only distinct values
        .collect()                               // Collect all distinct values
        .filter { sample -> sample.size() > 1 }  // Only proceed if more than 1 distinct value
        .flatMap { sample -> sample }            // Emit each distinct value separately
        .set{ ch_multisamples }
    // If there are multiple samples, make sure they have multiple assemblies
    ch_assembly
        .map { meta, file -> [meta.sample, meta, file] }
        .combine(ch_multisamples)                                                          // When ch_multisamples is empty, nothing will be forwarded
        .filter { sample_name, meta, files, valid_sample -> sample_name == valid_sample }  // The previous step produced too many combinations, reduce to genuine entries
        .map { sample_name, metas, files, valid_sample -> [sample_name, metas, files] }
        .groupTuple(by: 0) // Group by samples
        .map { sample_name, _meta, files -> [ [id: sample_name], files ] } // Drop meta information
        .filter { meta, files -> files.size() > 1 } // Only keep samples that have several assemblies
        .set { ch_to_quast_bysample }

    if(params.skip_kmerfinder){
        QUAST(
            ch_to_quast,
            params.reference_fasta ? [[:], reference_fasta] : [[:],[]],
            params.reference_gff ? [[:], reference_gff] : [[:],[]]
        )
        ch_quast_multiqc = QUAST.out.results
        QUAST_BYSAMPLE(
            ch_to_quast_bysample,
            params.reference_fasta ? [[:], reference_fasta] : [[:],[]],
            params.reference_gff ? [[:], reference_gff] : [[:],[]]
        )
    } else if (!params.skip_kmerfinder) {
        // Quast runs twice if kmerfinder is allowed.
        QUAST(
            ch_to_quast,
            [[:],[]],
            [[:],[]]
        )
        QUAST_BYSAMPLE(
            ch_to_quast_bysample,
            [[:],[]],
            [[:],[]]
        )
        QUAST_BYREFSEQID(
            ch_to_quast_byrefseq.map{ refmeta, consensus, _ref_fasta, _ref_gff -> tuple( refmeta, consensus)},
            ch_to_quast_byrefseq.map{ refmeta, _consensus, ref_fasta, _ref_gff -> tuple( refmeta, ref_fasta)},
            ch_to_quast_byrefseq.map{ refmeta, _consensus, _ref_fasta, ref_gff -> tuple( refmeta, ref_gff)}
        )
        ch_quast_multiqc = QUAST_BYREFSEQID.out.results
    }

    // Check assemblies that require further processing for gene annotation
    ch_assembly
        .branch{ _meta, fasta ->
            gzip: fasta.name.endsWith('.gz')
            skip: true
        }
        .set{ ch_assembly_for_gunzip }

    //
    // MODULE: BUSCO, assess genome assembly completeness
    //
    ch_busco_multiqc = channel.empty()
    if (!params.skip_busco) {
        BUSCO_BUSCO (
            ch_assembly,                                                        // tuple val(meta), path(fasta)
            params.busco_mode,                                                  // val mode
            params.busco_lineage,                                               // val lineage
            params.busco_db_path ? file(params.busco_db_path) : [],             // path busco_lineages_path
            params.busco_config_file ? file(params.busco_config_file) : [],     // path config_file (optional)
            params.busco_clean_intermediates                                    // val clean_intermediates
        )
        ch_busco_multiqc = BUSCO_BUSCO.out.short_summaries_txt
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions)
    }

    //
    // MODULE: PROKKA, gene annotation
    //
    ch_prokka_txt_multiqc = channel.empty()
    if ( !params.skip_annotation && params.annotation_tool == 'prokka' ) {
        // Uncompress assembly for annotation if necessary
        GUNZIP ( ch_assembly_for_gunzip.gzip )
        ch_to_prokka    = ch_assembly_for_gunzip.skip.mix( GUNZIP.out.gunzip )
        ch_versions     = ch_versions.mix( GUNZIP.out.versions )

        PROKKA (
            ch_to_prokka.filter{ _meta, fasta -> !fasta.isEmpty() },
            ch_proteins,
            []
        )
        ch_prokka_txt_multiqc   = PROKKA.out.txt.map{ _meta, prokka_txt -> [ prokka_txt ]}
        ch_versions             = ch_versions.mix(PROKKA.out.versions)
    }

    //
    // MODULE: BAKTA, gene annotation
    //
    ch_bakta_txt_multiqc = channel.empty()
    if ( !params.skip_annotation && params.annotation_tool == 'bakta' ) {
        // Uncompress assembly for annotation if necessary
        GUNZIP ( ch_assembly_for_gunzip.gzip )
        ch_to_bakta     = ch_assembly_for_gunzip.skip.mix( GUNZIP.out.gunzip )
        ch_versions     = ch_versions.mix( GUNZIP.out.versions )

        BAKTA_DBDOWNLOAD_RUN (
            ch_to_bakta.filter{ _meta, fasta -> !fasta.isEmpty() },
            params.baktadb,
            params.baktadb_download
        )
        ch_bakta_txt_multiqc    = BAKTA_DBDOWNLOAD_RUN.out.bakta_txt_multiqc.map{ _meta, bakta_txt -> [ bakta_txt ]}
        ch_versions             = ch_versions.mix(BAKTA_DBDOWNLOAD_RUN.out.versions)
    }
    //
    // MODULE: DFAST, gene annotation
    //
    // TODO: "dfast_file_downloader.py --protein dfast --dbroot ." could be used in a separate process and the db could be forwarded
    if ( !params.skip_annotation && params.annotation_tool == 'dfast' ) {
        DFAST (
            ch_assembly,
            channel.value(params.dfast_config ? file(params.dfast_config) : "")
        )
        ch_versions = ch_versions.mix(DFAST.out.versions)
    }

    //
    // MODULE: LIFTOFF, protein annotation
    //
    if ( !params.skip_annotation && params.annotation_tool == 'liftoff' ) {
        if (params.skip_kmerfinder || !params.liftoff_ref_from_kmerfinder) {
            // check if the reference files (fasta, gff) are given
            if ( !params.reference_fasta || !params.reference_gff ) {
                log.error "ERROR: when using liftoff with user specified reference, the `params.reference_fasta` and `params.reference_gff` must be provided."
            }

            LIFTOFF (
                ch_assembly,
                reference_fasta,
                reference_gff,
                []
            )
            ch_versions = ch_versions.mix(LIFTOFF.out.versions)
        } else {
            // run liftoff with kmerfinder reference

            LIFTOFF (
                ch_to_quast_byrefseq.map{ refmeta, consensus, _ref_fasta, _ref_gff -> tuple( refmeta, consensus)},
                ch_to_quast_byrefseq.map{ refmeta, _consensus, ref_fasta, _ref_gff -> tuple( refmeta, ref_fasta)},
                ch_to_quast_byrefseq.map{ refmeta, _consensus, _ref_fasta, ref_gff -> tuple( refmeta, ref_gff)},
                []
            )
            ch_versions = ch_versions.mix(LIFTOFF.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'bacass_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = !params.skip_kmerfinder && params.assembly_type ? channel.fromPath("$projectDir/assets/multiqc_config_${params.assembly_type}.yml", checkIfExists: true) : channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? channel.fromPath(params.multiqc_methods_description, checkIfExists: true) : channel.fromPath("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    CUSTOM_MULTIQC (
        ch_multiqc_config.ifEmpty([]),
        ch_multiqc_custom_config.ifEmpty([]),
        ch_multiqc_logo.ifEmpty([]),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        ch_multiqc_custom_methods_description.ifEmpty([]),
        ch_collated_versions.ifEmpty([]),
        ch_fastqc_raw_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_fastqc_trim_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_fastp_json_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_nanoplot_txt_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_porechop_log_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_filtlong_log_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_pycoqc_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_kraken_short_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_kraken_long_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_busco_multiqc.collect{it -> it[1]}.ifEmpty([]),
        ch_prokka_txt_multiqc.collect().ifEmpty([]),
        ch_bakta_txt_multiqc.collect().ifEmpty([]),
        ch_kmerfinder_multiqc.collectFile(name: 'multiqc_kmerfinder.yaml').ifEmpty([]),
    )

    emit:
    multiqc_report = CUSTOM_MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                        // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
