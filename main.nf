#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { single_end_kneaddata; paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan;rename_metaphlan_database_version; metaphlan_bzip} from "${projectDir}/processes/metaphlan.nf"
include { humann; humann_rename } from "${projectDir}/processes/humann.nf"



    // main workflow
    workflow {
        
        if (params.paired_end == true){
            println "Running paired_end_workflow"
            read_ch = Channel.fromFilePairs("${params.readsdir}/${params.filepattern}")
            println "Running paired_end_workflow"
            knead_out     =         paired_end_kneaddata(read_ch)
            } 
            
        else if (params.paired_end == false){
            println "Running single_end_workflow"
            read_ch = Channel
                .fromPath("${params.readsdir}/${params.filepattern}")
                .map { file -> 
                    def sample = file.baseName  // ERR3405856.fastq -> ERR3405856
                    return tuple(sample, file)
                }

            knead_out     =         single_end_kneaddata(read_ch)

            }
            else {
            throw new Exception("The paired_end must be bool true or false, got '${params.paired_end}'")
            }
    
    metaphlan_out =         metaphlan(knead_out.sample, knead_out.kneads)
    rename_metaphlan_out =  rename_metaphlan_database_version(metaphlan_out.sample, metaphlan_out.profile, metaphlan_out.bowtie2, metaphlan_out.sam)
                            metaphlan_bzip(rename_metaphlan_out.sample, rename_metaphlan_out.sam_rename)
    humann_out =            humann(metaphlan_out.sample, knead_out.kneads, metaphlan_out.profile)
    humann_rename_out =     humann_rename(humann_out.sample, humann.genefamilies, humann.log, humann.reactions, humann.pathabundance)
    }
