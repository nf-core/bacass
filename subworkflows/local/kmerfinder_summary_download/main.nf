//
// Kmerfinder subworkflow for species identification & QC
//
include { UNTAR                            } from '../../../modules/nf-core/untar'
include { KMERFINDER_KMERFINDER            } from '../../../modules/local/kmerfinder/kmerfinder'
include { KMERFINDER_SUMMARY               } from '../../../modules/local/kmerfinder/summary'
include { KMERFINDER_FIND_WINNER_REFERENCE } from '../../../modules/local/kmerfinder/find_winner_reference'
include { NCBI_DATASETS_DOWNLOAD           } from '../../../modules/local/ncbi_datasets_download'

workflow KMERFINDER_SUMMARY_DOWNLOAD {
    take:
    reads                   // channel: [ meta, reads ]
    consensus               // channel: [ meta, consensus ]

    main:
    ch_versions = Channel.empty()

    // Prepare kmerfinder database
    ch_kmerfinderdb = file(params.kmerfinderdb, checkIfExists: true)

    if ( ch_kmerfinderdb.name.endsWith('.gz') ) {
        UNTAR ( [[ id: ch_kmerfinderdb.getSimpleName() ], ch_kmerfinderdb] )
        ch_kmerfinderdb_untar = UNTAR.out.untar.map{ _meta, file -> file }
        ch_versions = ch_versions.mix(UNTAR.out.versions)
    } else {
        ch_kmerfinderdb_untar = Channel.fromPath(ch_kmerfinderdb)
    }
    ch_kmerfinderdb_untar = ch_kmerfinderdb_untar.map { it -> it.toAbsolutePath() }

    KMERFINDER_KMERFINDER (
        reads,    // channel: [ meta, reads ]
        ch_kmerfinderdb_untar.collect(),
        'bacteria'           // Val: 'tax_group'
    )
    ch_kmerfinder_report    = KMERFINDER_KMERFINDER.out.report
    ch_kmerfinder_json      = KMERFINDER_KMERFINDER.out.json
    ch_versions             = ch_versions.mix(KMERFINDER_KMERFINDER.out.versions)

    // MODULE: Kmerfinder summary report. Generates a csv report file collecting all sample references.
    KMERFINDER_SUMMARY (
        ch_kmerfinder_report.map{ _meta, report -> report }.collect()
    )
    ch_summary_yaml     = KMERFINDER_SUMMARY.out.yaml
    ch_versions         = ch_versions.mix(KMERFINDER_SUMMARY.out.versions)

    // SUBWORKFLOW: Create a channel to organize assemblies and reports based on the identified Kmerfinder reference.
    ch_kmerfinder_json
        .join(ch_kmerfinder_report, by:0)
        .join(consensus, by:0)
        .map{
            meta, report_json, report_txt, fasta ->
                def species_hits = report_json.splitJson(path:"kmerfinder.results.species_hits").value
                def specie = species_hits.size() > 0 ? species_hits.get(0)["Species"] : "Unknown Species"

                return tuple(specie, meta, report_txt, fasta)
        }
        .groupTuple(by:0) // Group by the "Species" field
        .set { ch_reports_byreference }

    // MODULE: Find the winner reference for each species
    KMERFINDER_FIND_WINNER_REFERENCE (
        ch_reports_byreference
            .map{ specie, _meta, report_txt, _fasta -> tuple(specie, report_txt) }
            .filter{ specie, _report_txt -> specie != "Unknown Species" }
    )
    ch_versions = ch_versions.mix(KMERFINDER_FIND_WINNER_REFERENCE.out.versions)

    // Prepare channel for NCBI_DATASETS_DOWNLOAD
    // Extract base accession from winner file (remove assembly version)
    ch_accessions_for_download = KMERFINDER_FIND_WINNER_REFERENCE.out.winner
        .map { refmeta, winner_file ->
            def full_accession = winner_file.text.trim()
            // Extract base accession: GCF_002795805.1_ASM279580v1 → GCF_002795805.1
            def base_accession = full_accession.split('_')[0] + '_' + full_accession.split('_')[1]
            return tuple([id: base_accession], base_accession)
        }

    // MODULE: Download reference genomes using NCBI datasets CLI
    NCBI_DATASETS_DOWNLOAD (
        ch_accessions_for_download
    )
    ch_versions = ch_versions.mix(NCBI_DATASETS_DOWNLOAD.out.versions)

    // Organize sample assemblies into channels based on their corresponding reference files.
    ch_reports_byreference
        .map { specie, meta, report_txt, fasta ->
            // Extract base accession from the first report to match with downloads
            def first_line = report_txt[0].text.split('\n').find { !it.startsWith('#') && it.trim() }
            def full_accession = first_line ? first_line.split('\t')[0] : null
            def base_accession = full_accession ? full_accession.split('_')[0] + '_' + full_accession.split('_')[1] : null
            return tuple(base_accession, specie, meta, report_txt, fasta)
        }
        .filter { base_accession, specie, meta, report_txt, fasta -> base_accession != null }
        .join(
            NCBI_DATASETS_DOWNLOAD.out.fna.map { meta, fna -> tuple(meta.id, fna) },
            by: 0
        )
        .join(
            NCBI_DATASETS_DOWNLOAD.out.gff.map { meta, gff -> tuple(meta.id, gff) },
            by: 0
        )
        .map {
            base_accession, specie, meta, report_txt, fasta, fna, gff ->
                return tuple([id: base_accession], meta, fasta, fna, gff)
        }
        .set { ch_consensus_byrefseq }

    emit:
    versions            = ch_versions               // channel: [ path(versions.yml) ]
    summary_yaml        = ch_summary_yaml           // channel: [ path(kmerfinder_summary.yml) ]
    consensus_byrefseq  = ch_consensus_byrefseq     // channel: [ refmeta, meta, fasta, fna, gff ]
}
