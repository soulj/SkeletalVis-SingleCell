#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_DOWNLOAD } from '../../modules/bam_download'

workflow test_BAM_DOWNLOAD {

    sampleTable = file("tests/testData/test_sampleTable.txt")
    ref = BAM_DOWNLOAD(sampleTable)
}
