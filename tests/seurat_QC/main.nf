#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { SEURAT_QC } from '../../modules/seurat_QC'




workflow test_SEURAT_QC {

  def input = []
    input = ["testSample1",file("tests/testData/testCounts")]

    t2g = file("tests/testData/transcripts_to_genes.txt")

    metaData = file("tests/testData/metaData.txt")

    SEURAT_QC(input,t2g,metaData)
}
