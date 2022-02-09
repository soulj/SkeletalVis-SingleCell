process SEURAT_MARKERGENES{
    tag "seurat marker genes $sampleName"
    label 'big_mem'


    publishDir "${params.outdir}/seurat_${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path (seurat) 
    
    output:
    tuple val(sampleName), path ("*_seurat_cluster.RDS")  , emit: seurat
    tuple val(sampleName), path ("*_seurat_cluster_slim.RDS") 
    tuple val(sampleName), path ("*.txt")   , emit: markers
    path ("*.png")

    script:
    """
    Rscript ${params.scriptDir}/seurat_markergenes.R \
    --seuratObject=$seurat \
    --seuratOutput=${sampleName}_seurat_cluster.RDS \
    --seuratDietOutput=${sampleName}_seurat_cluster_slim.RDS \
    --markerPlot=${sampleName}_markers.png \
    --diffExpTable="${sampleName}_markergenes.txt"
    """

}
