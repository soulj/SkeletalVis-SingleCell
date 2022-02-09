#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { SEURAT_MARKERGENES } from '../../modules/seurat_markergenes'

workflow test_SEURAT_MARKERGENES {

  def input = []
    input = ["testSample1",file("tests/testData/testSampleA_seruat_cellCycle.RDS")]
 
    SEURAT_MARKERGENES(input)
}
