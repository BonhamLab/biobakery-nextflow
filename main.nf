#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { single_end_kneaddata, paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; rename_metaphlan_database_version; metaphlan_bam } from "${projectDir}/processes/metaphlan.nf"
include { humann; humann_regroup; humann_rename } from "${projectDir}/processes/humann.nf"

println "readsdir: ${params.readsdir}"
println "filepattern: ${params.filepattern}"

read_ch.view { tuple -> "Running kneaddata on sample ${tuple[0]}, file ${tuple[1]}" }

// this workflow takes read_ch (a tuple of sample names and file names) as input 
workflow {

    if (params.paired_end == True){
        read_ch = Channel
        .fromFilePairs("${params.readsdir}/${params.filepattern}")
        } else if (params.paired_end == False){
        read_ch = Channel
            .fromPath("${params.readsdir}/${params.filepattern}")
            .map { file -> 
                def sample = file.baseName  // ERR3405856.fastq -> ERR3405856
                return tuple(sample, file)
            }
        }
        else {
        throw new Exception("The paired_end must be bool True or False, got '${params.paired_end}'")
        }
    knead_out     = kneaddata(read_ch)
    metaphlan_out = metaphlan(knead_out.sample, knead_out.fastq)
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
}


