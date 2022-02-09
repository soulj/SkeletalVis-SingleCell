process SRA_DOWNLOAD {

    label 'multi_cpu'

    publishDir path: "${params.outdir}/fastqFiles", mode: 'copy'

    input:
    tuple val (sampleID), val(SRR), val(Condition)

    output:
    path ('*R?.fastq.gz'), emit: fastqFiles

    script:
    
    """
    #set the sra cache
    vdb-config -s "/repository/user/main/public/root=$baseDir"

    #download the sra file and validate
    prefetch $SRR
    vdb-validate $SRR

    #dump the sra file into the three fastq files
    parallel-fastq-dump --split-files --gzip --sra-id $SRR --threads $params.cpuCores \
    --tmpdir ./

    #change the file names to R1,R2,I1 expectations if needed
    if [ -f ${SRR}_1.fastq.gz ]; then
       mv ${SRR}_1.fastq.gz  ${SRR}_I1.fastq.gz 
    fi

     if [ -f ${SRR}_2.fastq.gz ]; then
       mv ${SRR}_2.fastq.gz  ${SRR}_R1.fastq.gz 
    fi

     if [ -f ${SRR}_3.fastq.gz ]; then
       mv ${SRR}_3.fastq.gz  ${SRR}_R2.fastq.gz 
    fi
    """

}