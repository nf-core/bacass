//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    // reconfigure channels
    reads
        .map { meta, reads, long_fastq, fast5 -> [ meta, reads ] }
        .filter{ meta, reads -> reads != 'NA' } // TODO: filter when PE and 'NA'
        .set { shortreads }
    reads
        .map { meta, reads, long_fastq, fast5 -> [ meta, long_fastq, fast5 ] }
        .filter{ meta, long_fastq, fast5 -> long_fastq != 'NA' && fast5 != 'NA' }
        .set { longreads }

    emit:
    reads      // channel: [ val(meta), [ reads ], long_fastq, fast5 ]
    shortreads // channel: [ val(meta), [ reads ] ]
    longreads  // channel: [ val(meta), long_fastq, fast5 ]  
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ], long_fastq, fast5 ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.genome_size  = row.genome_size == null ? 'NA' : row.genome_size

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    // check long_fastq
    if (!row.long_fastq == 'NA') {
        if (!file(row.long_fastq).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Long FastQ file does not exist!\n${row.fastq_1}"
        }
        long_fastq = file(row.long_fastq)
    } else { long_fastq = 'NA' }

    // check long_fastq
    if (!row.fast5 == 'NA') {
        if (!file(row.fast5).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Fast5 file does not exist!\n${row.fastq_1}"
        }
        fast5 = file(row.fast5)
    } else { fast5 = 'NA' }

    // prepare output
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ], long_fastq, fast5 ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], long_fastq, fast5 ]
    }
    return array
}
