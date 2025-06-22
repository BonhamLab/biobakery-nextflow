# Nextflow pipeline for running the bioBakery

by Kevin Bonham, PhD 

bioBakery

- `KneadData`: a data quality-control pipeline that removes host genomic data within our metagenomic samples. Particularly, this pipeline uses a database containing a reference human genome so that all human DNA is removed from the samples. Link to more information here: (https://huttenhower.sph.harvard.edu/kneaddata/).
- `MetaPhlAn`
- `HUMAnN`

## Setup
Instructions for setting up a local environment to run the pipeline can be found on Danielle's notebook [here](LINK TO BE ADDED). 

Computing environments on the Tufts HPC and AWS should already be set-up with apptainer environments.

## Running the pipeline
This nextflow pipeline can be run on three different types of machines: 
1) Locally
2) Tufts high performance cluster (HPC)
3) Amazon website services cloud (AWS)

Based on the profiles described in `nextflow.config`, we can run the pipeline with the following Nextflow commands:

[NEED TO DOUBLE CHECK THIS PART]
### Running locally
`nextflow main.nf --local` 

### Running on the HPC
TO DO: Still need to figure out the exact nextflow syntax

Jobs on the Tufts HPC can be run in two different ways:
- **Batch**: the job will be sent to the queue and it will be completed based on how many resources you have requested, current cluster load, and fairshare (have you recently used the cluster) 
    - `nextflow main.nf --tufts_hpc --batch` 

- **Preempt**: this allows you to run your job preemptively using free nodes from another lab that paid for these compute resources. However, if they are already running a job, your job will be killed and you'll have to resubmit it.

    - `nextflow main.nf --tufts_hpc --preempt` 


### Running on AWS
`nextflow main.nf --amazon` 

> Kevin may want to add additional comments here about different ways to run the pipeline

## Databases
Several databases must be installed to run the pipeline. 

### Kneaddata
- A database containing a reference human genome so that unwanted human DNA can be removed from our metagenomic samples.

### MetaPhlAn
- `mpa_vOct22_CHOCOPhlAnSGB_202403` is the most recent MetaPhlAn database that is compatible with the versions of HUMAnN we are using
- Note: there is a more up-to-date version (released in January 2025) that we will probably eventually want to shift to once HUMAnN is able to support it.

### HUMAnN


## Information on software versions
This pipeline supports the following versions of MetaPhlAn and HUMAnN:
 
 ### MetaPhlAn
- MetaPhlAn 3.1.0
- MetaPhlAn 4

### HUMAnN
- HUMAnN3 3.7
- HUMAnN3 4 alpha

## Testing the pipeline
There are some raw fastq files in `test/` which can be processed through the pipeline

## Using the `master-params.yaml` file
The `master-params.yaml` file defines all input parameters that you may want to use to run the Nextflow pipeline. The file should not be used directly to run the pipeline. Rather, the user should select the params they need from the file based on how they would like to use the pipeline (software versions of MetaPhlAn or HUMAnN, computing environment, databases, etc. ), and paste these into a separate yaml file. This second yaml file can be used to run the Nextflow pipeline. 

### Overview of parameters in `master-params.yaml`
- `paired_end`: True or False, given the type of input data
- `metaphlan_ver`: MetaPhlAn software version (either `metaphlan3.1.0` or `metaphlan4`)
- `humann_ver`: HUMAnN3 software version (either `humann3.7` or `humann4_alpha`)
- `readsdir`: path to directory that contains raw data (bam files)
- `outdir`: path to directory where processed results will be saved
- `human_genome`: path to directory that contains human reference database used during Kneaddata 
- `metaphlan_db`: path to directory that contains metaphlan databases
- `metaphlan_index`: 
- `humann_nucleotide_db`: 
- `humann_protein_db`: 
- `humann_utility_db`: 
- `filepattern`: regex describing samples should be named (relative to the input raw data)