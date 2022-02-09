#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { SEURAT_CELLCYCLENORMALISE } from '../../modules/seurat_cellCycleNormalise'

workflow test_SEURAT_CELLCYCLENORMALISE {

  def input = []
    input = ["testSample1",file("tests/testData/testSampleA_seruat_qc.RDS")]
 
    SEURAT_CELLCYCLENORMALISE(input)
}
