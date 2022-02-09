process BAM_DOWNLOAD {

  label 'multi_cpu'

  input:
  path sampleTable

  output:
  path ('*.bam'), emit: bams

  shell:

  '''
  #get the File column from the sampleTable
  for file in `tail -n +2 !{sampleTable} | cut -d$'\t' -f 2`
  do

    #get the metadata xml
    efetch -db sra -id $file -format xml > metadata.xml

    #parse the xml to extract the url and md5 for the bam file
    xmlstarlet sel -t -m "//SRAFiles/SRAFile" -v "@url" \
    -o "checksum=md5=" -v "@md5" -n metadata.xml | \
    grep "bam" >> urls.txt

    #add the output file name
    echo " out=${file}.bam" >> urls.txt

  done

  #format the options for aria2c
  sed -i "s/checksum/\\n checksum/g" urls.txt
 

  #download all the bam files
  aria2c --input-file urls.txt
  '''

}
