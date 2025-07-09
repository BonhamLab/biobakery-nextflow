#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kneaddata } from './processes/kneaddata.nf'
include { metaphlan; rename_metaphlan_database_version } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

// if there is only 1 fastq per sample (not paired-end data)
workflow {
    println "readsdir: ${params.readsdir}"
    println "filepattern: ${params.filepattern}"
    read_ch = Channel
        .fromPath("${params.readsdir}/${params.filepattern}")
        .map { file -> 
            def sample = file.baseName  // ERR3405856.bam -> ERR3405856
            return tuple(sample, file)
        }

    read_ch.view { tuple -> "Running kneaddata on sample ${tuple[0]}, file ${tuple[1]}" }

    knead_out     = kneaddata(read_ch)
    metaphlan_out = metaphlan(knead_out.sample, knead_out.fastq)
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
}


