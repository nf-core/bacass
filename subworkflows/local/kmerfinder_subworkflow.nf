//
// Kmerfinder subworkflow for species identification & QC
//
include { UNTAR                     } from '../../modules/nf-core/untar/main'
include { KMERFINDER                } from '../../modules/local/kmerfinder/main'
include { KMERFINDER_SUMMARY        } from '../../modules/local/kmerfinder_summary/main'
include { FIND_DOWNLOAD_REFERENCE   } from '../../modules/local/find_download_reference'

workflow KMERFINDER_SUBWORKFLOW {
    take:
    reads                   // channel: [ meta, reads ]
    consensus               // channel: [ meta, consensus ]

    main:
    ch_versions = Channel.empty()

    // Prepare kmerfinder database
    ch_kmerfinderdb           = file(params.kmerfinderdb, checkIfExists: true)
    ch_ncbi_assembly_metadata = file(params.ncbi_assembly_metadata, checkIfExists: true)

    if ( ch_kmerfinderdb.name.endsWith('.gz') ) {
        UNTAR ( [[ id: ch_kmerfinderdb.getSimpleName() ], ch_kmerfinderdb] )
        ch_kmerfinderdb_untar = UNTAR.out.untar.map{ meta, file -> file }
        ch_versions = ch_versions.mix(UNTAR.out.versions)
    } else {
        ch_kmerfinderdb_untar = Channel.from(params.kmerfinderdb)
    }

    // MODULE: Kmerfinder, QC for sample purity. Identifies reference specie and reference genome assembly for each sample.
    reads
        .combine(ch_kmerfinderdb_untar)
        .map{ meta, reads, db -> tuple(meta, reads, db) }
        .set{ ch_to_kmerfinder }

    KMERFINDER (
        ch_to_kmerfinder,    // Channel: [ meta, reads, path_to_kmerfinderdb ]
        'bacteria'           // Val: 'tax_group'
    )
    ch_kmerfinder_report    = KMERFINDER.out.report
    ch_kmerfinder_json      = KMERFINDER.out.json
    ch_versions             = ch_versions.mix(KMERFINDER.out.versions)

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
                specie = report_json.splitJson(path:"kmerfinder.results.species_hits").value.get(0)["Species"]
                return tuple(specie, meta, report_txt, fasta)
        }
        .groupTuple(by:0) // Group by the "Species" field
        .set { ch_reports_byreference }

    // SUBWORKFLOW: For each species target, this subworkflow collects reference genome assemblies ('GCF*') and subsequently downloads the best matching reference assembly.
    FIND_DOWNLOAD_REFERENCE (
        ch_reports_byreference.map{ specie, meta, report_txt, fasta-> tuple(specie, report_txt) },
        ch_ncbi_assembly_metadata
    )
    ch_versions = ch_versions.mix(FIND_DOWNLOAD_REFERENCE.out.versions)

    // Organize sample assemblies into channels based on their corresponding reference files.
    ch_reports_byreference
        .join(FIND_DOWNLOAD_REFERENCE.out.fna)
        .join(FIND_DOWNLOAD_REFERENCE.out.gff)
        .join(FIND_DOWNLOAD_REFERENCE.out.winner)
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
