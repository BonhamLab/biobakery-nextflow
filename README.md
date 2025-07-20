# Nextflow pipeline for running the bioBakery

by Kevin Bonham, PhD and Danielle Pinto

[bioBakery](https://github.com/biobakery): software, documentation, and tutorials for microbial community profiling (created and mantained by the Huttenhower lab)

- [`KneadData`](https://github.com/biobakery/kneaddata): 
  a data quality-control pipeline that trims low quality reads
  and removes host genomic data within our metagenomic samples.
  Particularly, this pipeline uses a database containing a reference human genome
  so that all human DNA is removed from the samples.
  Link to more information here: (https://huttenhower.sph.harvard.edu/kneaddata/).
- [`MetaPhlAn`](https://github.com/biobakery/MetaPhlAn): 
  a computational tool for species-level microbial profiling (bacteria, archaea, eukaryotes, and viruses)
  from metagenomic shotgun sequencing data.
  Link to more information here:(https://huttenhower.sph.harvard.edu/metaphlan)
- [`HUMAnN`](https://github.com/biobakery/humann): 
  a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways
  in a community from metagenomic or metatranscriptomic sequencing data
  (typically millions of short DNA/RNA reads).
  This process, referred to as functional profiling,
  aims to describe the metabolic potential of a microbial community and its members.
  Link to more information here:(https://huttenhower.sph.harvard.edu/humann)

## Environment setup
Instructions for setting up a local environment to run the pipeline can be found on Danielle's notebook [here](https://github.com/BonhamLab/daniellepinto/blob/main/PeriodicMeetings/2025-06-17.md#danielles-personal-notes). 

Computing environments on the Tufts HPC and AWS should already be set-up with container-based (docker, apptainer) or conda environments.

## Running the pipeline
This nextflow pipeline can be run on three different types of machines: 
1) Locally
2) Tufts high performance cluster (HPC)
3) Amazon website services cloud (AWS)

Based on the profiles described in `nextflow.config`, we can run the pipeline with the following Nextflow commands:


### Running locally
`nextflow run main.nf --profile local -params-file template-params.yaml` 

> Note: To be able to run the pipeline, you must have at least 15 GB of available memory (RAM). This is needed to run the memory-intensive mapping step in metaphlan. Anecdotally, a laptop with 16 GB RAM was insufficient to run the pipeline.



### Running on the HPC

Jobs on the Tufts HPC can be run in two different ways:

- **Batch**: the job will be sent to the queue
  and it will be completed based on how many resources you have requested,
  current cluster load,
  and fairshare (have you recently used the cluster) 

- **Preempt**: this allows you to run your job using free nodes from another lab that paid for these compute resources.
  However, if they attempt to queue a job, your job will be preempted and killed, so you'll have to resubmit it.

With how the HPC environment is currently defined in `nextflow.config`,
jobs will first be submitted to the `batch` or `preempt` queue, whichever is available first.


- `nextflow run main.nf --profile tufts_hpc -params-file template-params.yaml` 

### Running on AWS
`nextflow main.nf --profile amazon -params-file template-params.yaml` 

> Note: We can also process samples on the MIT `engaging` cluster, but that should probably not be used without permission

## Databases
Several databases must be installed to run this pipeline. 

### Kneaddata
- A database containing a reference human genome so that unwanted human DNA can be removed from our metagenomic samples.
    - The `Homo_sapiens_hg39_T2T_Bowtie2_v0.1` bowtie2 database can be downloaded from [here](https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz).
        - This version of the database can be used for all analyses and there shouldn't be a big need to upgrade the database (unless we have an updated human genome!)
- Other reference databases can be added as well if other types of data want to be removed (eg. human transcriptome, mouse genome, etc.)

### MetaPhlAn
- `mpa_vOct22_CHOCOPhlAnSGB_202403` is the most recent MetaPhlAn database that is compatible with HUMAnN v4.
    - It can be found/downloaded manually from [here](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/). The easiest way to download is by running `metaphlan --install #any_other_args`
    - Note: there is a more up-to-date version (released in January 2025) that we will probably eventually want to shift to once HUMAnN is able to support it.
- `mpa_vOct22_CHOCOPhlAnSGB_202212` should be used if running HUMAnN v3.

#### A couple of pain points/things to know when downloading the Metaphlan databases:
- You cannot download multiple metaphlan databases to the same directory
- Downloading the bowtie2 database may take a bit of time (~30-60 minutes). If anything disrupts the download, remove the partially-downloaded file and try again. Any partially downloaded files will confuse metaphlan.
- MD5 files must be downloaded along with the bowtie databases in a single download. They should not be downloaded separately from the bowtie databases (which you may be tempted to do if your connection fails after downloading the bowtie databases), as the MD5 checksum must match the bowtie database checksum. 


### HUMAnN

- Use the command `humann_databases --download` to download the database that you desire 
    - for example: `humann_databases --download chocophlan full ~/biobakery_databases/chocophlan`
    - Each humann software version has its own matching database. As long as you are using the correct software version to run `humann_database --download`, you will get a compatible database. 
- All humann database can be found [here](http://cmprod1.cibio.unitn.it/databases/HUMAnN/).


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

## Using the `template-params.yaml` file
The `template-params.yaml` file defines all input parameters that you may want to use to run the Nextflow pipeline. The file should **not** be used directly to run the pipeline. Rather, the user should select the params they need from the file based on how they would like to use the pipeline (software versions of MetaPhlAn or HUMAnN, computing environment, databases, input data etc. ), and paste these into a separate yaml file. This second yaml file can be used to run the Nextflow pipeline. 

### Overview of parameters in `template-params.yaml`
- `paired_end`: True or False, given the type of input data
- `filepattern`: regex describing sample naming convention (relative to the input data type)
    - If there are paired-end reads, make sure the pattern considers both R1 and R2 the same sample.

- `metaphlan_version`: MetaPhlAn software version (either `metaphlan_v3` or `metaphlan_v4`)
- `humann_version`: HUMAnN3 software version (either `humann_v37` or `humann_v4a`)
- `readsdir`: path to directory that contains raw data 
- `outdir`: path to directory where processed results will be saved
- `human_genome`: path to directory that contains human reference database used during Kneaddata 
- `trimmomatic_path`: path to directory that contains Trimmomatic download
- `metaphlan_db`: path to directory that contains metaphlan databases
- `metaphlan_index`: database version (database must exist within `metaphlan_db`)
- `humann_db`: path to directory containing humann databases

