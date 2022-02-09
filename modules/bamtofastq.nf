process BAMTOFASTQ {

  label 'multi_bigmem'

  publishDir path: "${params.outdir}/fastqFiles", mode: 'copy'

  input:
  path bamFile

  output:
  path ('*R?.fastq.gz'), emit: fastqFiles

  script:

  """
  #run bamtofastq with reads-per-fastq set v.high to prevent fastq splitting
  bamtofastq --nthreads=${params.cpuCores} --reads-per-fastq=99999999999999999 $bamFile output

  #rename the fastq file to match the bam file
  mv output/*/*.fastq.gz .
  for x in *.fastq.gz; do
    mv "\$x" "${bamFile.baseName}_\${x%_*}.fastq.gz"
  done	
  """

}
