#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fastq } from './processes/bam2fastq.nf'
include { kneaddata } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_bam } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'


workflow {
    
    read_ch = Channel
        .fromFilePairs("${params.readsdir}/${params.filepattern}", size: 1)
	
 
    read_ch.view()
    
    bam_out       = bam2fastq(read_ch)
    knead_out     = kneaddata(bam_out)
    metaphlan_out = metaphlan(knead_out.sample, knead_out.fastq)
    metaphlan_bam = metaphlan_bam(metaphlan_out.sample, metaphlan_out.sam)
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
}
