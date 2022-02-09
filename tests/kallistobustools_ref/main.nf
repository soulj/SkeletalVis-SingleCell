#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.species = "human"
params.referenceDir = "$baseDir/../../../references"

include { KALLISTOBUSTOOLS_REF } from '../../modules/kallistobustools_ref'

workflow test_KALLISTOBUSTOOLS_REF {

    ref = KALLISTOBUSTOOLS_REF(params.species)
}
