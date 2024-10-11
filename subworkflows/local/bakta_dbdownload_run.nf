//
// Annotation of Bacterial genomes with Bakta
//

include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { BAKTA_BAKTA           } from '../../modules/nf-core/bakta/bakta/main'


workflow BAKTA_DBDOWNLOAD_RUN {
    take:
    ch_fasta                // channel: [ val(meta), path(fasta)  ]
    ch_path_baktadb         // channel: [ path(fasta) ]
    val_baktadb_download    // value: boolean

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Parse, download and/or untar Bakta database
    //
    if( ch_path_baktadb ){
        if (ch_path_baktadb.endsWith('.tar.gz')){
            ch_baktadb_tar  = Channel.from(ch_path_baktadb).map{ db -> [ [id: 'baktadb'], db ]}

            // MODULE: untar database
            UNTAR( ch_baktadb_tar )
            ch_baktadb      = UNTAR.out.untar.map{ meta, db -> db }
            ch_versions     = ch_versions.mix(UNTAR.out.versions)
        } else {
            ch_baktadb      = Channel.from(ch_path_baktadb).map{ db -> db }
        }
    } else if (!ch_path_baktadb && val_baktadb_download){
        // MODULE: Downlado Bakta database from zenodo
        BAKTA_BAKTADBDOWNLOAD()
        ch_baktadb  = BAKTA_BAKTADBDOWNLOAD.out.db
        ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)

    } else if (!ch_path_baktadb && !val_baktadb_download ){
        exit 1, "The Bakta database argument is missing. To enable the workflow to access the Bakta database, please include the path using '--baktadb' or use '--bakdtadb_download true' to download the Bakta database."
    }

    //
    // MODULE: BAKTA, gene annotation
    //
    BAKTA_BAKTA (
        ch_fasta,
        ch_baktadb,
        [],
        []
    )
    ch_bakta_txt_multiqc    = BAKTA_BAKTA.out.txt
    ch_versions             = ch_versions.mix(BAKTA_BAKTA.out.versions)

    emit:
    versions                = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
    bakta_txt_multiqc       = ch_bakta_txt_multiqc      // channel: [ meta, path(*.txt)  ]
}
