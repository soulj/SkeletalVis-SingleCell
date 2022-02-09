#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"
params.cellTypeDatabase = "HumanPrimaryCellAtlasData"


include { CELL_ANNOTATE } from '../../modules/cell_annotate'

workflow test_CELL_ANNOTATE {

  def input = []
    input = ["testSample1",file("tests/testData/testSampleA_seruat_cellCycle.RDS")]
 
    CELL_ANNOTATE(input)
}
