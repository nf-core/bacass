//
// Annotation of Bacterial genomes with Bakta
//

include { AUTOCYCLER_SUBSAMPLE }       from '../../../modules/nf-core/autocycler/subsample/main'
include { RAVEN as AUTOCYCLER_RAVEN }  from '../../../modules/nf-core/raven/main'
include { FLYE as AUTOCYCLER_FLYE }    from '../../../modules/nf-core/flye/main'
include { CANU as AUTOCYCLER_CANU }    from '../../../modules/nf-core/canu'
include { FASTA_CONSENSUS_AUTOCYCLER } from '../../../subworkflows/nf-core/fasta_consensus_autocycler/main'

workflow AUTOCYCLER {
    take:
    ch_preprocessed_fastq // channel: [ val(meta), path(short_reads), path(long_reads)  ]
    val_assemblers // channel: [ assembler ]
    val_flye_mode
    val_canu_mode

    main:
    ch_versions = channel.empty()
    ch_assemblies = channel.empty()

    // subsample and transpose to one subset per channel entry
    AUTOCYCLER_SUBSAMPLE ( 
        ch_preprocessed_fastq.map{ meta, _short_reads, long_reads -> [meta, long_reads] }, 
        ch_preprocessed_fastq.map { meta, _reads, _lr -> meta.gsize } 
    )
    AUTOCYCLER_SUBSAMPLE.out.subsampled_reads
        .transpose() // transpose to [ meta, fasta ]
        .map{ meta, reads -> 
            def new_meta = meta.clone()
            new_meta.subsample = reads.getBaseName() -'.fastq'
            [ new_meta, reads ]
        }
        .set{ ch_subsamples }

    // for each subset run chosen assemblers
    if ( val_assemblers.contains("raven") ) {
        AUTOCYCLER_RAVEN ( ch_subsamples )
        ch_assemblies = ch_assemblies.mix( AUTOCYCLER_RAVEN.out.fasta )
    }
    if ( val_assemblers.contains("flye") ) {
        AUTOCYCLER_FLYE (
            ch_subsamples,
            val_flye_mode
        )
        ch_assemblies = ch_assemblies.mix( AUTOCYCLER_FLYE.out.fasta )
    }
    if ( val_assemblers.contains("canu") ) {
        AUTOCYCLER_CANU (
            ch_subsamples,
            val_canu_mode,
            ch_subsamples.map { meta, _lr -> meta.gsize }
        )
        ch_assemblies = ch_assemblies.mix( AUTOCYCLER_CANU.out.assembly )
        ch_versions = ch_versions.mix(AUTOCYCLER_CANU.out.versions)
    }

    // unify assemblies with SUBWORKFLOW FASTA_CONSENSUS_AUTOCYCLER
    FASTA_CONSENSUS_AUTOCYCLER ( 
        ch_assemblies
            .map{ meta, assembly ->
                def new_meta = meta.clone()
                new_meta.remove("subsample")
                tuple( new_meta, assembly)
            }
            .filter{ meta, assembly -> assembly.countLines() > 1 } // keep only non-empty assembly files
            .groupTuple() // group to "[ val(meta), [ fasta, fasta, ... ] ]"
        )

    emit:
    versions                = ch_versions.ifEmpty([]) // channel: [ path(versions.yml) ]
    assembly                = FASTA_CONSENSUS_AUTOCYCLER.out.consensus_assembly // channel: [ val(meta), fasta ]
    assembly_graph          = FASTA_CONSENSUS_AUTOCYCLER.out.consensus_assembly_graph // channel: [ val(meta), gfa ]
}
