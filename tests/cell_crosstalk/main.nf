#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"
params.cellTypeDatabase = "HumanPrimaryCellAtlasData"
params.cellGroups = ["SingleR.labels","seurat_annotations"]

include { CELL_CROSSTALK  } from '../../modules/cell_crosstalk'

workflow test_CELL_CROSSTALK {

  def input = []
    input = ["testSample1",file("tests/testData/lianaTestData.RDS")]
 
    CELL_CROSSTALK(input,params.cellGroups)
}
