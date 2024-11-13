#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'

workflow {
    
    read_pairs_ch = Channel
        .fromFilePairs(
            [ "$params.readsdir/$params.filepattern",
              "$params.readsdir/*_kneaddata.fastq.gz" ],
            size:-1)

    metaphlan_db      = params.metaphlan_db

    metaphlan_out = metaphlan(read_pairs_ch, metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
}
