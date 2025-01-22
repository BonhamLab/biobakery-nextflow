#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { prefetch_and_split } from './processes/prefetch_and_split.nf'
include { kneaddata } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

workflow {
    
    // Step 1: Read the SRA accession list and create a channel
    sra_ch = Channel.fromPath(params.sra_list)
        .splitText()
        .map { it.trim() }

    // Step 2: Prefetch and split FASTQ files
    fastq_pairs_ch = prefetch_and_split(sra_ch)
    
    human_genome      = params.human_genome
    metaphlan_db      = params.metaphlan_db
    
    knead_out     = kneaddata(fastq_pairs_ch, human_genome)
    metaphlan_out = metaphlan(knead_out[0], metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
}
