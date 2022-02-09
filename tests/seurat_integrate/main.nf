#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { SEURAT_INTEGRATE } from '../../modules/seurat_integrate'

workflow test_SEURAT_INTEGRATE {

  def input = []
    input = [file("tests/testData/test1_seurat_integrate.RDS"),file("tests/testData/test2_seurat_integrate.RDS")]
 
    SEURAT_INTEGRATE(input)
}
