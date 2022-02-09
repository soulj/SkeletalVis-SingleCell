process MAKEREACTOMEPATHWAYS {
    tag "makeReactomePathways"
    label 'big_mem'

    publishDir path: "${params.referenceDir}/enrichment", mode: 'copy'
    
    input:
    val  species

    output:
    path "*.RDS" , emit: reactome

    script:
    """
    Rscript ${params.scriptDir}/makeReactomePathways.R \
    --species=$species \
    --reactomePathways="${params.species}_reactomePathways.RDS"
    """
    
    
}