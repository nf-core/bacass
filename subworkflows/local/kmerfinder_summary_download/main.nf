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
    ch_assembly             // channel: [ meta, assembly ]

    main:
    ch_versions = channel.empty()

    // Prepare kmerfinder database
    ch_kmerfinderdb = file(params.kmerfinderdb, checkIfExists: true)

    if ( ch_kmerfinderdb.name.endsWith('.gz') ) {
        UNTAR ( [[ id: ch_kmerfinderdb.getSimpleName() ], ch_kmerfinderdb] )
        ch_kmerfinderdb_untar = UNTAR.out.untar.map{ _meta, file -> file }
        ch_versions = ch_versions.mix(UNTAR.out.versions)
    } else {
        ch_kmerfinderdb_untar = channel.fromPath(ch_kmerfinderdb)
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
        .cross(ch_assembly) { it -> it[0].sample } // merge by meta.sample -> [[ meta, report_json, report_txt ],[ meta, assembly ]]
            .map { kmerfinder, assembly ->
                def meta              = assembly[0]
                def kmerfinder_json   = kmerfinder[1]
                def kmerfinder_report = kmerfinder[2]
                def assembly_fasta    = assembly[1]
                [ meta, kmerfinder_json, kmerfinder_report, assembly_fasta ] }
        .map{
            meta, report_json, report_txt, fasta ->
                def species_hits = report_json.splitJson(path:"kmerfinder.results.species_hits").value
                def species = species_hits.size() > 0 ? species_hits.get(0)["Species"] : "Unknown Species"
                [ species, meta, report_txt, fasta ]
        }
        .groupTuple(by:0) // Group by the "Species" field
        .set { ch_reports_byreference }

    // MODULE: Find the winner reference for each species
    KMERFINDER_FIND_WINNER_REFERENCE (
        ch_reports_byreference
            .map{ species, _meta, report_txt, _fasta ->
                [ species, report_txt ] }
            .filter{ species, _report_txt -> species != "Unknown Species" }
    )
    ch_versions = ch_versions.mix(KMERFINDER_FIND_WINNER_REFERENCE.out.versions)

    // Preserve the winner selected for each species. This key must also be used
    // when associating the downloaded reference with the sample assemblies.
    // Deriving an accession from an arbitrary report in the species group can
    // select a different reference and leave the downstream joins empty.
    ch_winner_references = KMERFINDER_FIND_WINNER_REFERENCE.out.winner
        .map { species, winner_file ->
            def full_accession = winner_file.text.trim()
            // Extract base accession: GCF_002795805.1_ASM279580v1 → GCF_002795805.1
            def base_accession = full_accession.split('_')[0] + '_' + full_accession.split('_')[1]
            return tuple(species, base_accession)
        }

    // Prepare channel for NCBI_DATASETS_DOWNLOAD.
    ch_accessions_for_download = ch_winner_references
        .map { _species, base_accession -> tuple([id: base_accession], base_accession) }

    // MODULE: Download reference genomes using NCBI datasets CLI
    NCBI_DATASETS_DOWNLOAD (
        ch_accessions_for_download
    )
    ch_versions = ch_versions.mix(NCBI_DATASETS_DOWNLOAD.out.versions)

    // Associate each species with the reference selected by
    // KMERFINDER_FIND_WINNER_REFERENCE and its downloaded files.
    ch_winner_references
        .map { species, base_accession -> tuple(base_accession, species) }
        .join(
            NCBI_DATASETS_DOWNLOAD.out.fna.map { meta, fna -> tuple(meta.id, fna) },
            by: 0
        )
        .join(
            NCBI_DATASETS_DOWNLOAD.out.gff.map { meta, gff -> tuple(meta.id, gff) },
            by: 0
        )
        .map {
            base_accession, species, fna, gff ->
                return tuple(species, [id: base_accession], fna, gff)
        }
        .set { ch_reference_by_species }

    // Add every assembly in a species group to that species' selected reference.
    ch_reports_byreference
        .join(ch_reference_by_species, by: 0)
        .map {
            _species, meta, _report_txt, fasta, refmeta, fna, gff ->
                return tuple(refmeta, meta, fasta, fna, gff)
        }
        .set { ch_consensus_byrefseq }

    emit:
    versions            = ch_versions               // channel: [ path(versions.yml) ]
    summary_yaml        = ch_summary_yaml           // channel: [ path(kmerfinder_summary.yml) ]
    consensus_byrefseq  = ch_consensus_byrefseq     // channel: [ refmeta, meta, fasta, fna, gff ]
}
