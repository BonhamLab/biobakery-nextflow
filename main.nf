#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kneaddata } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'


workflow {
    
    read_pairs_ch = Channel
        .fromFilePairs("$params.readsdir/$params.filepattern", size: 2)

    human_genome      = params.human_genome
    metaphlan_db      = params.metaphlan_db
    
    knead_out     = kneaddata(read_pairs_ch, human_genome)
    metaphlan_out = metaphlan(knead_out[0], knead_out[1], metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
}
