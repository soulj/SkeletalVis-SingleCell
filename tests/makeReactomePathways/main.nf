#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir = "data/${params.accessionNumber}"
params.referenceDir = "$baseDir/../../../references"
params.scriptDir = "../../../src"
params.species = "human"


include { MAKEREACTOMEPATHWAYS } from '../../modules/makeReactomePathways'

workflow test_MAKEREACTOMEPATHWAYS {

  MAKEREACTOMEPATHWAYS(params.species)

}
