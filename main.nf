#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fastq } from "${projectDir}/processes/bam2fastq.nf"
include { kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; metaphlan_bam } from "${projectDir}/processes/metaphlan.nf"
include { humann; humann_regroup; humann_rename } from "${projectDir}/processes/humann.nf"


workflow {
    
    read_ch = Channel
        .fromFilePairs("${params.readsdir}/${params.filepattern}")
    
    knead_out     = kneaddata(read_ch)
    metaphlan_out = metaphlan(knead_out.sample, knead_out.fastq)
    // metaphlan_bam = metaphlan_bam(metaphlan_out.sample, metaphlan_out.sam) // not working because of headers
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
}
