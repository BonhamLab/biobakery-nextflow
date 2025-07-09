#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { metaphlan } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

workflow {
    
    read_ch = Channel
        .fromPath("${params.readsdir}/*_R1.fastq.gz")
        .map { file ->
            def sample = file.baseName.replaceAll("_R1.fastq", "")
            def paired = [
                "${params.readsdir}/${sample}_R1.fastq.gz",
                "${params.readsdir}/${sample}_R2.fastq.gz"
            ]
            def unpaired = [
                "${params.readsdir}/${sample}_R1_unpaired.fastq.gz", 
                "${params.readsdir}/${sample}_R2_unpaired.fastq.gz"
            ]
            return [sample, paired, unpaired]
        }
    
    metaphlan_out = metaphlan(read_ch.map{it-> it[0]}, read_ch.map{it -> it[1]}, read_ch.map{it -> it[2]})
    humann_out    = humann(metaphlan_out.sample, metaphlan_out.concatenated, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
}
