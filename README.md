
# SkeletalVis-SingleCell

## Introduction

**SkeletalVis-SingleCell** is a bioinformatics pipeline for reproducible analyses of 10x Genomics single-cell RNA-sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a portable workflow tool to run tasks across multiple compute infrastructures. This pipeline uses a singularity container containing all the software needed to run the analysis, making installation simple and the results reproducible.

## Pipeline summary

The **SkeletalVis-SingleCell** pipeline takes a sample table and a parameter file defining the experiment as input. If not provided fastq files are automatically downloaded using the provided sample identifiers.

### Features:
(**a**) Download of fastq files either directly from ENA, via conversion of sra or bam files from SRA<br/>
(**b**)	Quantification using [`kallisto-bustools`](https://www.kallistobus.tools/) to produce cell x gene matrices<br/>
(**c**) Flexible filtering of [`empty droplets`](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), quality control and thresholding<br/>
(**d**) Normalisation and cell cycle effect removal<br/>
(**e**) Automatic cell type annotation with [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html)<br/>
(**f**) Clustering and visualisation with [`Seurat`](https://satijalab.org/seurat/)<br/>
(**g**) Marker gene identification and pathway analysis<br/>
(**h**) Cell crosstalk analysis of ligand-receptor predictions using [`liana`](https://github.com/saezlab/liana)<br/>
(**i**) Sample integration and differential expression analysis between conditions with [`miloR`](https://github.com/MarioniLab/miloR)<br/>

Analyses are run in parallel and in result of error you can resume with the `-resume` parameter to re-run the pipeline starting from the previous fault.

## Quick Start

### Analyse an example dataset

Try the pipeline on an example dataset (all inputs will be automatically downloaded): -

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

3. [`Configure`](https://www.nextflow.io/docs/latest/config.html) the resource profile for your HPC or local computer. A template for slurm schedulers is provided as an example in `nextflow.config`

4. Download the pipeline and test on the example dataset with a single command:

    ```console
     nextflow run skeletalvis.nf -profile slurm -params-file GSE152805.yaml -with-singularity library://jsoul/default/singlecell:latest
    ```
### Analyse your own data

1. Define the sampleTable

Create a tab seperated table with unique Sample names, SRR accession numbers (if download is needed) and any additional metadata e.g

|Sample|File|Condition|
| ---|---|---|
|Control_1|SRRXXX	|Control|
|Control_2|	SRRXXX	|Control|
|Treated_1|	SRRXXX	|Treated|
|Treated_2|	SRRXXX	|Treated|

2. Define the configuration

Most parameters are set to sensible defaults within the main nextflow script, with only 5 parameters required to be altered with typical use:

|Parameter|Description|Options|
| ---|---|---|
|accession|The GEO accession of the data - used to name output data and download fastq files||
|downloadSite|The site to download the raw data from if needed	|SRA, ENA, SRA_BAM|
|species|The species the reads originate from - used to create the kallisto bus index	|human, mouse|
|chemistry|The chemistry used for the 10x Genomics experiment	|10xv1, 10xv2, 10xv3|
|replciates|Does the experiment contain replicated treatments to perform differential expression analysis?|true, false|


Parameters should be defined within a yaml file. See `params/GSE152805.yaml` for an example.

3. Run the pipeline with your own parameters

    ```console
     nextflow run skeletalvis.nf -profile slurm -params-file ownData.yaml -with-singularity library://jsoul/default/skeletalvis-singlecell
    ```

### Testing modules
Modules can be tested using the [`pytest-workflow`](https://pypi.org/project/pytest-workflow/) framework. Module test directories within the `tests` folder contain a nextflow script and a configuration yaml file defining the test for each module.

1. Install pytest-workflow

    ```console
	conda install pytest-workflow
    ```

2. Run the tests - e.g to test the GSEA module

    ```console
	pytest --symlink --kwdof --tag gsea
    ```


