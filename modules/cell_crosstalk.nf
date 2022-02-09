process CELL_CROSSTALK{

    tag "cell crosstalk $sampleName"
    label 'big_mem'


    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sampleName), path (seurat)
    each groupBy
    
    output:
    path ("*.txt")
    

    script:
    """
    Rscript ${params.scriptDir}/liana.R \
    --seuratObject=$seurat \
    --groupBy=$groupBy \
    --interactions=${sampleName}_${groupBy}_interactions.txt
    """

}
