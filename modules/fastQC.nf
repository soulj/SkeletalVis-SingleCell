process FASTQC {

  tag "$sampleId"
  publishDir "${params.outdir}/fastqc", mode: 'copy',
  saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(sampleId), path(reads)

  output:
  path "*_fastqc.{zip,html}", emit: stats

  script:
  """
  fastqc -t 1 -q $reads
  """

}
