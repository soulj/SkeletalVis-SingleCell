process CELLRANGER_COUNT {

  tag "cellranger count"

  publishDir path: "${params.outdir}/cellranger", mode: 'copy'

  label 'multi_big_mem'


  input:

  tuple val(name), path(reads)
  path(transcriptome)

  output:

  tuple val(name), path ("${name}/outs/raw_feature_bc_matrix")  , emit: counts_unfiltered

  shell:
  '''
  cellranger count --id=!{name} \
                   --transcriptome=!{transcriptome} \
                   --fastqs=$(echo !{reads}|tr " " "\n"|sort|tr "\n" " ") \
                   --sample=!{name} \
                   --localcores=!{params.cpuCores} \
                   --localmem=!{params.mem}
  '''

}
