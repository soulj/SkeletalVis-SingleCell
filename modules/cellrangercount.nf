process CELLRANGER_COUNT {

  tag "cellranger count"

  publishDir path: "${params.outdir}/cellranger", mode: 'copy'

  label 'cellranger'


  input:

  tuple val(name), path(reads)
  path(transcriptome)

  output:
  path ("${name}/outs")
  tuple val(name), path ("${name}/outs/raw_feature_bc_matrix")  , emit: counts_unfiltered

  shell:
  '''
  mkdir !{name}_fastqFiles
  mv !{reads} !{name}_fastqFiles
  
  s="!{reads}"
  sampleName="${s%%_*}"
  
   cellranger count --id=!{name} \
                   --transcriptome=!{transcriptome} \
                   --fastqs=!{name}_fastqFiles \
                   --sample=$sampleName \
                   --localcores=!{task.cpus} \
                   --localmem=!{task.memory.giga} \
                   --no-bam
  '''

}
