process MULTIQC {

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    label 'big_mem'
    input:
    path ('fastqc/*')

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    multiqc .
    """

}

