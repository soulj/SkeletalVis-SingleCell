#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.cpuCores = 1
params.transcriptome = "tests/testData/cellranger_tiny_ref/"
params.mem = 10

include { CELLRANGER_COUNT } from '../../modules/cellrangercount.nf'

workflow test_CELLRANGER_COUNT {

  def input = []
    input = ["tinygex",[
    file("tests/testData/cellranger_tiny_fastq/tinygex_S1_L001_R1_001.fastq.gz"),
    file("tests/testData/cellranger_tiny_fastq/tinygex_S1_L001_R2_001.fastq.gz")]]

    transcriptome = file(params.transcriptome)

    CELLRANGER_COUNT(input,transcriptome)
}
