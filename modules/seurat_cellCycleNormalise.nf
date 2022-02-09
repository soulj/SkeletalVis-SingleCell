process SEURAT_CELLCYCLENORMALISE{
    tag "seurat cell cycle $sampleName"
    label 'big_mem'


    publishDir "${params.outdir}/seurat_${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path (seurat) 
    
    output:
    tuple val(sampleName), path ("*.RDS")  , emit: seurat
    path ("*.png")

    script:
    """
    Rscript ${params.scriptDir}/seurat_cellCycleNormalise.R \
    --seuratObject=$seurat \
    --seuratOutput=${sampleName}_seruat_cellCycle.RDS \
    --preNormPlot=${sampleName}_preCellCycleNorm.png \
    --postNormPlot=${sampleName}_postCellCycleNorm.png
    """

}
