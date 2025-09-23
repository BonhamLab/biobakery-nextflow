#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { single_end_kneaddata; paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; metaphlan_bzip} from "${projectDir}/processes/metaphlan.nf"
include { humann} from "${projectDir}/processes/humann.nf"



    // main workflow
    workflow {
        
        if (params.paired_end == true){
            println "Running paired_end_workflow"
            read_ch = Channel.fromFilePairs("${params.readsdir}/${params.filepattern}", checkIfExists: true)
            println "Running paired_end_workflow"
            knead_out     =         paired_end_kneaddata(read_ch)
            } 
            
        else if (params.paired_end == false){
            println "Running single_end_workflow"
            read_ch = Channel
                .fromPath("${params.readsdir}/${params.filepattern}")
                .map { file -> 
                    def sample = file.baseName.replaceFirst(/(\.fastq|\.fq)$/, '')  // ERR3405856.fastq -> ERR3405856
                    return tuple(sample, file)
                }

            knead_out     =         single_end_kneaddata(read_ch)

            }
            else {
            throw new Exception("The paired_end must be bool true or false, got '${params.paired_end}'")
            }
    
    metaphlan_out =         metaphlan(knead_out.sample, knead_out.kneads)
    metaphlan_bzip_out =    metaphlan_bzip(metaphlan_out.sample, metaphlan_out.sam)
    humann_out =            humann(metaphlan_out.sample, knead_out.kneads, metaphlan_out.profile)
    }
