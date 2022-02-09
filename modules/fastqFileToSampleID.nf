process FASTQTOSAMPLEID {

   label 'big_mem'

   input:
   tuple val(sampleID), path(samples)
   path sampleTable

   output:
   tuple stdout , path(samples)

   script:
   """
   python2 ${params.scriptDir}/fastqToSampleID.py $sampleTable $sampleID
   """
}