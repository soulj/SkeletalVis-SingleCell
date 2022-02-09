process KB_COUNT {

  tag "kallistobustools count"

  publishDir path: "${params.outdir}/kallisto", mode: 'copy'

  label 'multi_big_mem'


  input:

  tuple val(name), path(reads)
  path index
  path t2g
  val chemistry

  output:

  tuple val(name), path ("${name}_bus_output/counts_unfiltered")  , emit: counts_unfiltered
  path ("${name}_kallisto.log") , emit: kallisto_log

  shell:
  '''
  kb count -i !{index} -g !{t2g} \
  -x !{chemistry} \
  -t!{params.cpuCores} \
  -o !{name}_bus_output/ $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ") &>  !{name}_kallisto.log
  '''

}
