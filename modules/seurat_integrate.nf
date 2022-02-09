process SEURAT_INTEGRATE{
    tag "seurat integrate $sampleName"
    label 'big_mem'


    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path (seuratFiles) 
    
    output:
    path ("*.RDS") , emit: seurat
    path ("*.png")

    script:
    """
    Rscript ${params.scriptDir}/seurat_integrate.R \
    --seuratObjects="${seuratFiles}" \
    --seuratOutput=seurat_integrated.RDS \
    --tSNEPlot=integrated_tSNE.png 
    """

}
