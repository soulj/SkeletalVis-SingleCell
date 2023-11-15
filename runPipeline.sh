#!/bin/bash
# SLURM script to run cell ranger count

#SBATCH --mem=4G


module load singularity
module load adoptopenjdk
module load nextflow

accession=$1


echo "running cellranger pipeline for $accession"

nextflow run -profile slurm -w tmp/${accession}_tmp main.nf -params-file params/${accession}.yaml -with-singularity singlecell_latest.sif -entry CellRanger
