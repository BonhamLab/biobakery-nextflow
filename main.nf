#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fastq } from './processes/bam2fastq.nf'
include { kneaddata } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'


workflow {
    
    read_ch = Channel
        .fromPath("$params.readsdir/$params.filepattern")

    human_genome      = params.human_genome
    metaphlan_db      = params.metaphlan_db
    humann_bowtie_db  = params.humann_bowtie_db
    humann_protein_db = params.humann_protein_db
    humann_utility_db = params.humann_utility_db
    
    bam_out       = bam2fastq(read_ch)
    knead_out     = kneaddata(bam_out, human_genome)
    metaphlan_out = metaphlan(knead_out.fastq, metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out.sample, metaphlan_out.sam)
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile, humann_bowtie_db, humann_protein_db)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies, humann_utility_db)
    humann_rename(regroup_out, humann_utility_db)
}
