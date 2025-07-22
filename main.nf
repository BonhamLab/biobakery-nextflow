#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { single_end_kneaddata; paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; preprocess_paired_end_metaphlan; rename_metaphlan_database_version; metaphlan_bam } from "${projectDir}/processes/metaphlan.nf"
include { humann; humann_regroup; humann_rename } from "${projectDir}/processes/humann.nf"

println "readsdir: ${params.readsdir}"
println "filepattern: ${params.filepattern}"

// workflow for single-end data
workflow single_end_workflow{
    read_ch = Channel
            .fromPath("${params.readsdir}/${params.filepattern}")
            .map { file -> 
                def sample = file.baseName  // ERR3405856.fastq -> ERR3405856
                return tuple(sample, file)
            }

    knead_out     =         single_end_kneaddata(read_ch)

    metaphlan_out =         metaphlan(knead_out.sample, knead_out.fastq)
    rename_metaphlan_out =  rename_metaphlan_database_version(metaphlan_out.sample)
    metaphlan_bam_out =     metaphlan_bam(metaphlan_out.sample, metaphlan_out.bam)

    humann_out =            humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    humann_regroup_out =    humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename_out =     humann_rename(humann_regroup_out.sample, humann_regroup_out.ecs, 
                            humann_regroup_out.kos, humann_regroup_out.pfams)


}


// workflow for paired-end data

workflow paired_end_workflow{
    read_ch = Channel
        .fromFilePairs("${params.readsdir}/${params.filepattern}")

    knead_out     = paired_end_kneaddata(read_ch)
    metaphlan_out = metaphlan(knead_out.sample, knead_out.fastq)
    humann_out    = humann(metaphlan_out.sample, knead_out.fastq, metaphlan_out.profile)
    regroup_out   = humann_regroup(humann_out.sample, humann_out.genefamilies)
    humann_rename(regroup_out)
    



}


// main workflow
workflow {

    if (params.paired_end == true){
        paired_end_workflow
        } else if (params.paired_end == false){
        single_end_workflow
        }
        else {
        throw new Exception("The paired_end must be bool true or false, got '${params.paired_end}'")
        }

    read_ch.view { tuple -> "Running kneaddata on sample ${tuple[0]}, file ${tuple[1]}" }
}
