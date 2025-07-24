#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { single_end_kneaddata; paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; preprocess_paired_end_metaphlan; rename_metaphlan_database_version; metaphlan_bam } from "${projectDir}/processes/metaphlan.nf"
include { humann; humann_regroup; humann_rename } from "${projectDir}/processes/humann.nf"

// workflow for single-end data
workflow single_end_workflow{
    main:
        read_ch = Channel
                .fromPath("${params.readsdir}/${params.filepattern}")
                .map { file -> 
                    def sample = file.baseName  // ERR3405856.fastq -> ERR3405856
                    return tuple(sample, file)
                }

        read_ch.view { tuple -> "Running kneaddata on sample ${tuple[0]}, file ${tuple[1]}" }

        knead_out     =         single_end_kneaddata(read_ch)

        metaphlan_out =         metaphlan(knead_out.sample, knead_out.fastq)
        rename_metaphlan_out =  rename_metaphlan_database_version(metaphlan_out.sample)
        metaphlan_bam_out =     metaphlan_bam(metaphlan_out.sample, metaphlan_out.bam)

        humann_out =            humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
        humann_regroup_out =    humann_regroup(humann_out.sample, humann_out.genefamilies)
        humann_rename_out =     humann_rename(humann_regroup_out.sample, humann_regroup_out.ecs, 
                                humann_regroup_out.kos, humann_regroup_out.pfams)

    emit:
    knead_out
}


// workflow for paired-end data

workflow paired_end_workflow {
  main: 
  {
    println "💥 paired_end_workflow was called"
  
    read_ch = Channel.fromFilePairs("${params.readsdir}/${params.filepattern}", flat: true)
    read_ch.view { tuple -> println "Running kneaddata on sample ${tuple[0]}, file ${tuple[1]}" }
   
    knead_out     =         paired_end_kneaddata(read_ch)
    metaphlan_out =         metaphlan(knead_out.sample, knead_out.fastq, knead_out.paired, knead_out.unpaired)
    rename_metaphlan_out =  rename_metaphlan_database_version(metaphlan_out.sample)
    metaphlan_bam_out =     metaphlan_bam(metaphlan_out.sample, metaphlan_out.bam)

    humann_out =            humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    humann_regroup_out =    humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename_out =     humann_rename(humann_regroup_out.sample, humann_regroup_out.ecs,
                            humann_regroup_out.kos, humann_regroup_out.pfams)
}}


    // main workflow
    workflow {
        println "🧪 params.paired_end = ${params.paired_end} (${params.paired_end.getClass()})"
        
        if (params.paired_end == true){
            println "Running paired_end_workflow"
            result = paired_end_workflow

            } else if (params.paired_end == false){
            println "Running single_end_workflow"
            result = single_end_workflow

            }
            else {
            throw new Exception("The paired_end must be bool true or false, got '${params.paired_end}'")
            }
    }
