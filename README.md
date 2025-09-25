# Nextflow pipeline for running the bioBakery

![Static Badge](https://img.shields.io/badge/Author-Kevin_Bonham_PhD-purple)
![Static Badge](https://img.shields.io/badge/Author-Guilherme_Fahur_Bottino_PhD-purple)
![Static Badge](https://img.shields.io/badge/Author-Danielle_Pinto-purple)

[bioBakery](https://github.com/biobakery): software, documentation,
and tutorials for microbial community profiling (created and mantained by the Huttenhower lab)

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

> [!Note]
> To be able to run the pipeline, you must have at least 15 GB of available memory (RAM). This is needed to run the memory-intensive mapping step in metaphlan. Anecdotally, a laptop with 16 GB RAM was insufficient to run the pipeline.
> Additionally, it has been notoriously challenging to download the metaphlan databases. It is recommended to run this pipeline on the Tufts HPC, where the environment and databases have already been established. 

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

- Prior to scheduling the job, make sure to load nextflow `module load nextflow` and add bin to your path (`export PATH="/cluster/tufts/bonhamlab/shared/bin:$PATH"`)

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
  database that is compatible with HUMAnN v4.
    - It can be found/downloaded manually from [here](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/).
      The easiest way to download is by running `metaphlan --install #any_other_args`
    - Note: there is a more up-to-date version (released in January 2025)
      that we will probably eventually want to shift to once HUMAnN is able to support it.
- `mpa_vOct22_CHOCOPhlAnSGB_202212` should be used if running HUMAnN v3.

#### A couple of pain points/things to know when downloading the Metaphlan databases:

- You cannot download multiple metaphlan databases to the same directory
- Downloading the bowtie2 database may take a bit of time (~30-60 minutes). If anything disrupts the download, remove the partially-downloaded file and try again. Any partially downloaded files will confuse metaphlan.
- MD5 files must be downloaded along with the bowtie databases in a single download. They should not be downloaded separately from the bowtie databases (which you may be tempted to do if your connection fails after downloading the bowtie databases), as the MD5 checksum must match the bowtie database checksum. 
- Due to these issues, it is recommended to run the pipeline on the Tufts HPC (where database files have already been downloaded), rather than trying to re-download them locally.

### HUMAnN

- Use the command `humann_databases --download` to download the database that you desire 
    - for example: `humann_databases --download chocophlan full ~/biobakery_databases/chocophlan`
    - Each humann software version has its own matching database. As long as you are using the correct software version to run `humann_database --download`, you will get a compatible database. 
- All humann database can be found [here](http://cmprod1.cibio.unitn.it/databases/HUMAnN/).

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

## Pipeline outputs

> [!NOTE] Aspirational
> Not all of this section may be reflected in the current pipeline

The general structure of outputs is

- `${UID}` is the unique sequence identifier, eg `SEQ99999_S1`
- `${DB}` is a chocophlan database, eg `mpa_vOct22_CHOCOPhlAnSGB_202403`
- `${H_VERSION}` is a `humann` version,
  eg `v4a` or `v37`

```txt
./raw_data_root_folder/
├── rawfastq/
    ├── ${UID}_l001_r1_0001.fastq.gz
    └── ${UID}_l001_r2_0001.fastq.gz

./processing_root_folder/
    ├── kneaddata/
    │   ├── ${UID}_kneaddata.log
    │   ├── ${UID}_kneaddata_paired_1.fastq.gz
    │   ├── ${UID}_kneaddata_paired_2.fastq.gz
    │   ├── ${UID}_kneaddata_unmatched_1.fastq.gz
    │   └── ${UID}_kneaddata_unmatched_2.fastq.gz
    ├── metaphlan
    │   └── ${DB}
    │       ├── ${UID}_${DB}_profile.tsv
    │       ├── ${UID}_${DB}_bowtie2.tsv
    │       └── ${UID}_${DB}.sam.bz2
    └── humann
        └── ${H_VERSION}
            ├── ${UID}_${H_VERSION}_genefamilies.tsv
            ├── ${UID}_${H_VERSION}_pathabundance.tsv
            └── ${UID}_${H_VERSION}_pathcoverage.tsv
```

> [!NOTE] Humann file formats
> In humann v4a, the outputs are slightly different.
> The files are structured like:
> 
> ```
>     └── humann
>        └── ${H_VERSION}
>            ├── ${UID}_${H_VERSION}_0.log
>            ├── ${UID}_${H_VERSION}_2_genefamilies.tsv
>            ├── ${UID}_${H_VERSION}_3_reactions.tsv
>            └── ${UID}_${H_VERSION}_4_pathabundance.tsv
> ```

### Deprecated structure

- `humann` used to have 3 subfolders, `main`, `regroup`, `rename`.
  These have been deprecated in favor of only saving primary outputs
  in long-term storage. 
  Regrouped and renamed files can be trivially created from primary outputs.


  ### Testing

  Test have been written using [`NF-test`](https://www.nf-test.com/) to test the end-to-end pipeline, as well as the outputs. To run the tests on the HPC cluster, you can run:
  - `nf-test test tests/main.nf.test --profile tufts_hpc`
  - `nf-test test tests/validate_output.nf.test --profile tufts_hpc`
  
  
   CI has also been set up on Github to ensure that new changes to the pipeline are tested automatically. 



