#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

workflow {
    // Channel for input files: read pre-existing kneaddata output files from a specified directory
    knead_ch = Channel
        .fromPath("$params.kneaddata_dir/*_kneaddata*.fastq.gz")
        .map { file ->
            def sample = file.name.replaceAll(/_kneaddata(?:_[paired|unmatched]_[1|2])?\.fastq\.gz$/, "")
            tuple(sample, file)
        }
        .groupBy { it[0] }  // Group by sample name to keep related files together

    // Define databases
    metaphlan_db      = params.metaphlan_db
    humann_bowtie_db  = params.humann_bowtie_db
    humann_protein_db = params.humann_protein_db
    humann_utility_db = params.humann_utility_db
    
    // Run metaphlan process directly on pre-existing kneaddata files
    metaphlan_out = metaphlan(knead_ch, metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[3])

}