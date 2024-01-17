//
// Kmerfinder subworkflow for species identification & QC
//
include { KMERFINDER                } from '../../modules/local/kmerfinder'
include { KMERFINDER_SUMMARY        } from '../../modules/local/kmerfinder_summary'
include { FIND_DOWNLOAD_REFERENCE   } from '../../modules/local/find_download_reference'
include { QUAST                     } from '../../modules/nf-core/quast/main'

workflow KMERFINDER_SUBWORKFLOW {
    take:
    kmerfinder_db           // channel: [ path ]
    ncbi_assembly_metadata  // channel: [ path ]
    reads                   // channel: [ meta, reads ]
    consensus               // channel: [ meta, consensus ]

    main:
    ch_versions = Channel.empty()

    // MODULE: Kmerfinder, QC for sample purity. Identifies reference specie and reference genome assembly for each sample.
    KMERFINDER (
        reads,
        kmerfinder_db
    )
    ch_kmerfinder_report    = KMERFINDER.out.report
    ch_kmerfinder_json      = KMERFINDER.out.json
    ch_versions             = ch_versions.mix( KMERFINDER.out.versions.ifEmpty(null) )

    // MODULE: Kmerfinder summary report. Generates a csv summary file collecting all sample reports.
    KMERFINDER_SUMMARY (
        ch_kmerfinder_report.map{meta, report -> report }.collect()
    )
    ch_summary_yaml     = KMERFINDER_SUMMARY.out.yaml
    ch_versions         = ch_versions.mix( KMERFINDER_SUMMARY.out.versions.ifEmpty(null) )

    // SUBWORKFLOW: Grouping reports by identified reference species.
    ch_kmerfinder_json
        .join(ch_kmerfinder_report, by:0)
        .join(consensus, by:0)
        .map{
            meta, report_json, report_txt, fasta ->
                specie = report_json.splitJson(path:"kmerfinder.results.species_hits").value.get(0)["Species"]
                return tuple(specie, meta, report_txt, fasta)
        }
        .groupTuple(by:0) // Group by the "Species" field
        .set { ch_reports_byreference }

    // SUBWORKFLOW: For each specie target, this subworkflow collects the reference genomes assemblies ('GCF*'), and subsequently, downloads the winning reference assembly.
    if (!params.reference_fasta && !params.reference_gff) {
        FIND_DOWNLOAD_REFERENCE (
            ch_reports_byreference.map{ specie, meta, report_txt, fasta-> tuple(specie, report_txt) },
            ncbi_assembly_metadata
        )
        ch_versions         = ch_versions.mix( FIND_DOWNLOAD_REFERENCE.out.versions.ifEmpty(null) )

        // Arrange sample's assemblyes into channels with their corresponding reference files.
        ch_reports_byreference
            .join(FIND_DOWNLOAD_REFERENCE.out.fna)
            .join(FIND_DOWNLOAD_REFERENCE.out.gff)
            .join(FIND_DOWNLOAD_REFERENCE.out.winner)
            .map {
                specie, meta, report_txt, fasta, fna, gff, winner_id ->
                    return tuple([id: winner_id.getBaseName()], meta, fasta, fna, gff)
            }
            .set { ch_consensus_byrefseq }

    } else if (params.reference_fasta && params.reference_gff) {
        // TODO: Haven't tested so far
        //ch_reports_byreference
        //    .join(params.reference_fasta)
        //    .join(params.reference_gff)
        //    .map {
        //        refmeta, meta, report_txt, fasta ->
        //            refmeta.id = params.reference_fasta.getBaseName()
        //            return tuple(refmeta, meta, fasta, fna, gff)
        //    }
        //    .set { ch_consensus_byrefseq  }
    }

    emit:
    versions            = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
    summary_yaml        = ch_summary_yaml           // channel: [ path(kmerfinder_summary.yml) ]
    consensus_byrefseq  = ch_consensus_byrefseq     // channel: [ refmeta, meta, fasta, fna, gff ]
}
