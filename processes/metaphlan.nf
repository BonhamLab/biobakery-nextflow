process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "*.tsv"

    input:
    tuple val(sample), path(kneads)
    path unmatched
    path metaphlan_db

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_bowtie2.tsv" , emit: bowtie2
    path "${sample}.sam"         , emit: sam


    script:
    """
    metaphlan $kneads ${sample}_profile.tsv \
        --mapout ${sample}_bowtie2.tsv \
        --samout ${sample}.sam \
        --input_type fastq \
        --nproc ${task.cpus} \
        --bowtie2db $metaphlan_db
    """
}
 
 process metaphlan_bam {
    tag "metaphlan_bam on $sample"
    publishDir "$params.outdir/metaphlan"
    stageInMode "copy"

    input:
    val sample
    path sam

    output:
    val  sample          , emit: sample
    path "${sample}.bam" , emit: bam

    when:

    script:
    """
    samtools -bS $sam -o ${sample}.bam
    """
}