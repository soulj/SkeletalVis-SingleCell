process GSEA{

    tag "gsea $sampleName"
    label 'big_mem'


    publishDir "${params.outdir}/seurat_${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path (markers)
    path(reactomePathways)
    
    output:
    path ("*.txt")


    script:
    """
    Rscript ${params.scriptDir}/gsea.R --diffExp=$markers \
    --reactomePathways=$reactomePathways \
    --enrichmentTable=${sampleName}_enrichmentTable.txt

    """

}