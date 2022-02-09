process SEURAT_QC{
    tag "seurat qc $sampleName"
    label 'big_mem'

    
    publishDir "${params.outdir}/seurat_${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path (countsDir) 
    path t2g
    path metaData
    
    output:
    tuple val(sampleName), path ("*.RDS")  , emit: seurat
    path ("*.png")

    script:
    """
    Rscript ${params.scriptDir}/seurat_qc.R \
    --tenx=$countsDir --t2g=$t2g \
    --metaData=$metaData \
    --sampleName="${sampleName}" \
    --seuratOutput=${sampleName}_seruat_qc.RDS \
    --QCPlot=${sampleName}_qc.png
    """
    
}
