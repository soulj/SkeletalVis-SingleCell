process MILO {

    tag "milo"
    label 'big_mem'


    publishDir "${params.outdir}/MILO", mode: 'copy'

    input:
    path (seu)
    
    output:
    path ("*.txt")
    path ("*.png")


    script:
    """
    Rscript ${params.scriptDir}/milo.R --seuratObject=$seu \
    --upRegulated=upRegulatedMarkers.txt \
    --downRegulated=downRegulatedMarkers.txt \
    --neighbourhoods=neighbourhoods.png

    """

}