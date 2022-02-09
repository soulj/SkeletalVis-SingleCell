#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../src"


include { MILO } from '../../modules/milo'

workflow test_MILO {

  input = file("tests/testData/test_integratedSeurat.RDS")
  
  MILO(input)
}
