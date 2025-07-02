//
// Kmerfinder subworkflow for species identification & QC
//
include { UNTAR                         } from '../../../modules/nf-core/untar'
include { KMERFINDER_KMERFINDER         } from '../../../modules/local/kmerfinder/kmerfinder'
include { KMERFINDER_SUMMARY            } from '../../../modules/local/kmerfinder/summary'
include { KMERFINDER_DOWNLOAD_REFERENCE } from '../../../modules/local/kmerfinder/download_reference'

workflow KMERFINDER_SUMMARY_DOWNLOAD {
    take:
    reads                   // channel: [ meta, reads ]
    consensus               // channel: [ meta, consensus ]

    main:
    ch_versions = Channel.empty()

    // Prepare kmerfinder database
    ch_kmerfinderdb           = file(params.kmerfinderdb, checkIfExists: true)

    if ( ch_kmerfinderdb.name.endsWith('.gz') ) {
        UNTAR ( [[ id: ch_kmerfinderdb.getSimpleName() ], ch_kmerfinderdb] )
        ch_kmerfinderdb_untar = UNTAR.out.untar.map{ meta, file -> file }

        ch_versions = ch_versions.mix(UNTAR.out.versions)
    } else {
        ch_kmerfinderdb_untar = Channel.fromPath(ch_kmerfinderdb)
    }
    ch_kmerfinderdb_untar = ch_kmerfinderdb_untar.map { it -> it.toAbsolutePath() }

    KMERFINDER_KMERFINDER (
        reads,    // Channel: [ meta, reads ]
        ch_kmerfinderdb_untar.collect(),
        'bacteria'           // Val: 'tax_group'
    )
    ch_kmerfinder_report    = KMERFINDER_KMERFINDER.out.report
    ch_kmerfinder_json      = KMERFINDER_KMERFINDER.out.json
    ch_versions             = ch_versions.mix(KMERFINDER_KMERFINDER.out.versions)

    // MODULE: Kmerfinder summary report. Generates a csv report file collecting all sample references.
    KMERFINDER_SUMMARY (
        ch_kmerfinder_report.map{ meta, report -> report }.collect()
    )
    ch_summary_yaml     = KMERFINDER_SUMMARY.out.yaml
    ch_versions         = ch_versions.mix(KMERFINDER_SUMMARY.out.versions)

    // SUBWORKFLOW:  Create a channel to organize assemblies and reports based on the identified Kmerfinder reference.
    ch_kmerfinder_json
        .join(ch_kmerfinder_report, by:0)
        .join(consensus, by:0)
        .map{
            meta, report_json, report_txt, fasta ->
                species_hits = report_json.splitJson(path:"kmerfinder.results.species_hits").value
                def specie = species_hits.size() > 0 ? species_hits.get(0)["Species"] : "Unknown Species"

                return tuple(specie, meta, report_txt, fasta)
        }
        .groupTuple(by:0) // Group by the "Species" field
        .set { ch_reports_byreference }

    // SUBWORKFLOW: For each species target, this subworkflow collects reference genome assemblies ('GCF*') and subsequently downloads the best matching reference assembly.
    KMERFINDER_DOWNLOAD_REFERENCE (
        ch_reports_byreference
            .map{ specie, meta, report_txt, fasta-> tuple(specie, report_txt) }
            .filter{ specie, report_txt -> specie != "Unknown Species" }
    )
    ch_versions = ch_versions.mix(KMERFINDER_DOWNLOAD_REFERENCE.out.versions)

    // Organize sample assemblies into channels based on their corresponding reference files.
    ch_reports_byreference
        .join(KMERFINDER_DOWNLOAD_REFERENCE.out.fna)
        .join(KMERFINDER_DOWNLOAD_REFERENCE.out.gff)
        .join(KMERFINDER_DOWNLOAD_REFERENCE.out.winner)
        .map {
            specie, meta, report_txt, fasta, fna, gff, winner_id ->
                return tuple([id: winner_id.getBaseName()], meta, fasta, fna, gff)
        }
        .set { ch_consensus_byrefseq }

    emit:
    versions            = ch_versions               // channel: [ path(versions.yml) ]
    summary_yaml        = ch_summary_yaml           // channel: [ path(kmerfinder_summary.yml) ]
    consensus_byrefseq  = ch_consensus_byrefseq     // channel: [ refmeta, meta, fasta, fna, gff ]
}
