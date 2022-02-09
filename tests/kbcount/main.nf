#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.chemistry="10xv2"
params.cpuCores = 1

include { KB_COUNT } from '../../modules/kbcount'

workflow test_KB_COUNT {

  def input = []
    input = ["testSample1",
    [file("https://github.com/pachterlab/kb_python/raw/master/tests/fixtures/R1.fastq"),file("https://github.com/pachterlab/kb_python/raw/master/tests/fixtures/R2.fastq")]]

    t2g=file("https://github.com/pachterlab/kb_python/raw/master/tests/fixtures/transcripts_to_genes.txt")

     index= file("https://github.com/pachterlab/kb_python/raw/master/tests/fixtures/mouse_truncated.idx")

    KB_COUNT(input,index,t2g,params.chemistry)
}
