//
// Kmerfinder subworkflow for species identification & QC
//
include { KMERFINDER                } from '../../modules/local/kmerfinder'
include { KMERFINDER_SUMMARY        } from '../../modules/local/kmerfinder_summary'
include { FIND_DOWNLOAD_REFERENCE   } from '../../modules/local/find_download_reference'
include { QUAST                     } from '../../modules/nf-core/quast/main'

workflow KMERFINDER_SUBWORKFLOW {
    take:
    kmerfinder_db       // channel: [ path ]
    ncbi_assembly_metadata    // channel: [ path ]
    reads               // channel: [ meta, reads ]
    consensus           // channel: [ meta, consensus ]

    main:
    ch_versions = Channel.empty()

    // MODULE: Kmerfinder, QC for sample purity
    KMERFINDER (
        reads,
        kmerfinder_db
    )
    ch_kmerfinder_report    = KMERFINDER.out.report
    ch_kmerfinder_json      = KMERFINDER.out.json
    ch_versions             = ch_versions.mix( KMERFINDER.out.versions.ifEmpty(null) )

    // MODULE: Kmerfinder summary report
    KMERFINDER_SUMMARY (
        ch_kmerfinder_report.map{meta, report -> report }.collect()
    )
    ch_summary_yaml     = KMERFINDER_SUMMARY.out.yaml
    ch_versions         = ch_versions.mix( KMERFINDER_SUMMARY.out.versions.ifEmpty(null) )

    // SUBWORKFLOW: Group assemblies by reference geneome
    ch_kmerfinder_json
        .join(ch_kmerfinder_report, by:0)
        .join(consensus, by:0)
        .map{
            meta, report_json, report_txt, fasta ->
                def refseq = [:]
                refseq.id = report_json.splitJson(path:"kmerfinder.results.species_hits").value.get(0)["Assembly"]
                return tuple(refseq, meta, report_txt, fasta)
        }
        .groupTuple(by:0)
        .set { ch_consensus_byrefseq }

    // MODULE: Find & Download common reference sequences
    if (!params.reference_fasta && !params.reference_gff) {
        FIND_DOWNLOAD_REFERENCE (
            ch_consensus_byrefseq.map{ refseq, meta, report_txt, fasta -> tuple(refseq, report_txt)},
            ncbi_assembly_metadata
        )
        ch_reference_fasta  = FIND_DOWNLOAD_REFERENCE.out.fna
        ch_reference_gff    = FIND_DOWNLOAD_REFERENCE.out.gff
        ch_versions         = ch_versions.mix( FIND_DOWNLOAD_REFERENCE.out.versions.ifEmpty(null) )
    }

    // Get reference sequence IDs
    ch_consensus_byrefseq
        .map{ refseq, meta, report_txt, fasta -> refseq }
        .collect()
        .set { ch_refseqid }

    emit:
    versions            = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
    summary_yaml        = ch_summary_yaml           // channel: [ path(kmerfinder_summary.yml) ]
    refseqids           = ch_refseqid               // channel: [ val(refseq1), val(refseq1),...]
    reference_fasta     = ch_reference_fasta        // channel: [ meta,  path(*.fna) ]
    reference_gff       = ch_reference_gff          // channel: [ meta,  path(*.gff) ]
    consensus_byrefseq  = ch_consensus_byrefseq     // channel: [ refseq, meta, report_txt, fasta ]
}
