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

Instructions for setting up a local environment to run the pipeline
can be found on Danielle's notebook [here](https://github.com/BonhamLab/daniellepinto/blob/main/PeriodicMeetings/2025-06-17.md#danielles-personal-notes). 

Computing environments on the Tufts HPC and AWS should already be set-up with container-based
(docker, apptainer) or conda environments.

## Running the pipeline

This nextflow pipeline can be run on three different types of machines: 
1) Locally
2) Tufts high performance cluster (HPC)
3) Amazon website services cloud (AWS)

Based on the profiles described in `nextflow.config`,
we can run the pipeline with the following Nextflow commands:


### Running locally

```sh
nextflow run main.nf -profile local -params-file params.yaml`
```

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

```sh
nextflow run main.nf -profile tufts_hpc -params-file params.yaml
```

### Running on AWS

The AWS setup is too complicated to describe here for the moment,
but requires setting up and correctly configuring:

- Amazon machine images (AMIs) with the correct software installed and resources available
  - note that often the IOPS (I/O operations per second) can be overwhelmed if multiple
    jobs are run on the same image and all try to pull containers at the same time
- Job queues that define which AMI should be used for different jobs
- machine users with the correct permissions
- S3 buckets with correct permissions
- Containers with the correct versions of different software uploaded to ECR or docker-hub

```sh
nextflow main.nf -profile amazon -params-file params.yaml`
```

## Databases

Several databases must be installed to run this pipeline. 

### Kneaddata

- A database containing a reference human genome so that unwanted human DNA can be removed from our metagenomic samples.
    - The `Homo_sapiens_hg39_T2T_Bowtie2_v0.1` bowtie2 database
      can be downloaded from [here](https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz).
        - This version of the database can be used for all analyses
          and there shouldn't be a big need to upgrade the database (unless we have an updated human genome!)
- Other reference databases can be added as well if other types
  of data want to be removed (eg. human transcriptome, mouse genome, etc.)

### MetaPhlAn

- `mpa_vOct22_CHOCOPhlAnSGB_202403` is the most recent MetaPhlAn
  database that is compatible with the versions of HUMAnN we are using.
    - It can be found/downloaded manually from [here](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/).
      The easiest way to download is by running `metaphlan --install #any_other_args`
- Note: there is a more up-to-date version (released in January 2025) that we will probably
  eventually want to shift to once HUMAnN is able to support it.

### HUMAnN

The `humann_databases` command should be used to interact with and download
the databases used by the software.
Minimally there are 3 databases required:

- `nucleotide` - this defines the link between ChocoPhlAn marker genes database
  and genomes of taxa that `humann` will first attempt to align to
- `protein` - this is a DIAMOND (by default) database to align protein sequences
  if the nucleotide-based search does not find a hit
- `utility_mapping` - This defines things like mappings between UniRef90 and other annotations (KO, EC, etc),
  human-readable names for database accessions, etc.

The available databases can be seen with using the following command:

```sh
humann_databases --available
```

To download a specific database,
use `humann_databases --download` and include the database and the particular build, as well as the download location
For example:

```sh
humann_databases --download chocophlan full /some/path/to/databases
```

By default, this will also update the local configuration.
You can use `--upadate-config no` do disable this.

An alternate approach is to download databases manually,
eg [http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz].
To set the `humann` configuration to point to an existing download,
you can use `humann_config`, eg:

```sh
humann_config --update database_folders nucleotide /some/path/to/databases/chocophlan
```

The location of databases can also be set independently of the configuration
by passing command line flags when running `humann`.
For example `humann --nucleotide-database /some/path/to/databases/chocophlan ...`

## Information on software versions

This pipeline supports the following versions of MetaPhlAn and HUMAnN:
 
 ### MetaPhlAn

- MetaPhlAn v3.1.0
- MetaPhlAn v4

### HUMAnN

- HUMAnN v3.7
- HUMAnN v4 alpha

## Testing the pipeline

There are some raw fastq files in `test/` which can be processed through the pipeline

## Using the `template-params.yaml` file

The `template-params.yaml` file defines all input parameters that you may want to use to run the Nextflow pipeline. The file should **not** be used directly to run the pipeline. Rather, the user should select the params they need from the file based on how they would like to use the pipeline (software versions of MetaPhlAn or HUMAnN, computing environment, databases, input data etc. ), and paste these into a separate yaml file. This second yaml file can be used to run the Nextflow pipeline. 

### Overview of parameters in `template-params.yaml`

- `input_data_type`: type of input data (either `fastq` or `bam`)
- `paired_end`: True or False, given the type of input data
- `filepattern`: regex describing sample naming convention (relative to the input data type)
- `metaphlan_version`: MetaPhlAn software version (either `metaphlan_v3` or `metaphlan_v4`)
- `humann_version`: HUMAnN3 software version (either `humann_v37` or `humann_v4a`)
- `readsdir`: path to directory that contains raw data 
- `outdir`: path to directory where processed results will be saved
- `human_genome`: path to directory that contains human reference database used during Kneaddata 
- `metaphlan_db`: path to directory that contains metaphlan databases
- `metaphlan_index`: database version (database must exist within `metaphlan_db`)
- `humann_nucleotide_db`: path to directory containing chocophlan database
- `humann_protein_db`: path to directory containing UniRef database
- `humann_utility_db`: path to directory containing databases that have conversions
  between different protein annotations (eg UniRef90 to KO or EC),
  and names for all of the different annotations that have them

