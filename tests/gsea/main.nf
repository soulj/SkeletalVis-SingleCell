#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir = "data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { GSEA } from '../../modules/gsea'

workflow test_GSEA {

  def input = []
    input = ["testSample1",file("tests/testData/testSample1_markergenes.txt")]

  reactomePathways = file("tests/testData/human_reactomePathways.RDS")


  GSEA(input,reactomePathways)

}
