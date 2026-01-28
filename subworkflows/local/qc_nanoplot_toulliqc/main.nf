/*
 * Nanopore Long Read QC with Nanoplot and ToulligQC
 */

include { NANOPLOT  } from '../../../modules/nf-core/nanoplot/main'
include { TOULLIGQC } from '../../../modules/nf-core/toulligqc/main'

workflow QC_NANOPLOT_TOULLIGQC {

    take:
    ch_fastq
    skip_nanplot
    skip_toulligqc

    main:

    /*
     * Nanopore QC with Nanoplot
     */

    nanoplot_png     = channel.empty()
    nanoplot_html    = channel.empty()
    nanoplot_txt     = channel.empty()
    nanoplot_log     = channel.empty()
    nanoplot_version = channel.empty()
    if (!skip_nanplot) {
        NANOPLOT ( ch_fastq )
        nanoplot_png     = NANOPLOT.out.png
        nanoplot_html    = NANOPLOT.out.html
        nanoplot_txt     = NANOPLOT.out.txt
        nanoplot_log     = NANOPLOT.out.log
        nanoplot_version = NANOPLOT.out.versions
    }

    /*
     * Nanopore QC with ToulligQC
     */
    toulligqc_report_data   = channel.empty()
    toulligqc_report_html   = channel.empty()
    toulligqc_plots_html    = channel.empty()
    toulligqc_plotly_js     = channel.empty()
    toulligqc_version       = channel.empty()
    if (!skip_toulligqc) {
        TOULLIGQC ( ch_fastq )
        toulligqc_report_data  = TOULLIGQC.out.report_data
        toulligqc_report_html  = TOULLIGQC.out.report_html
        toulligqc_plots_html   = TOULLIGQC.out.plots_html
        toulligqc_plotly_js    = TOULLIGQC.out.plotly_js
        toulligqc_version      = TOULLIGQC.out.versions
    }

    emit:
    nanoplot_png
    nanoplot_html
    nanoplot_txt
    nanoplot_log
    nanoplot_version

    toulligqc_report_data
    toulligqc_report_html
    toulligqc_plots_html
    toulligqc_plotly_js
    toulligqc_version
}
