process CELL_ANNOTATE{
    tag "cell_annotate $sampleName"
    label 'big_mem'

    
    publishDir "${params.outdir}/seurat_${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path (seurat) 
    
    output:
    tuple val(sampleName), path ("*.RDS")  , emit: seurat
    path ("*.png")
    path ("*.txt")

    script:
    """
    Rscript ${params.scriptDir}/singleR_annotate.R \
    --seuratObject=$seurat \
    --seuratOutput=${sampleName}_seurat_annotated.RDS \
    --cellTypeDatabase=${params.cellTypeDatabase} \
    --cellTypePlot=${sampleName}_cellAnnotations.png \
    --cellTypeStats=${sampleName}_cellStats.txt \
    --cellTypeMarkers=${sampleName}_cellTypeMarkers.txt
    """
    
}
