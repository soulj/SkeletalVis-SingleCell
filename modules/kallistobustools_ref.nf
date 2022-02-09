process KALLISTOBUSTOOLS_REF {
    tag "kallistobustools ref"
    label 'big_mem'

    publishDir path: "${params.referenceDir}", mode: 'copy'
    
    input:
    val  species

    output:
    path "*_Index.idx" , emit: index
    path "*_t2g.txt"   , emit: t2g


    script:
    """
    kb ref -d $species -g ${params.species}_t2g.txt -i ${params.species}_Index.idx
    """
    
    
}
