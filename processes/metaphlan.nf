process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "*.tsv"

    input:
    val(sample)
    path(kneads)

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_bowtie2.tsv" , emit: bowtie2
    path "${sample}.sam"         , emit: sam


    script:
    """
    metaphlan $kneads -o ${sample}_profile.tsv \
        --mapout ${sample}_bowtie2.tsv \
        --samout ${sample}.sam \
        --input_type fastq \
        --nproc ${task.cpus} \
        --db_dir ${params.metaphlan_db} \
        --index ${params.metaphlan_index} \
        -t rel_ab_w_read_stats
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
    path "${sample}_markers.bam" , emit: bam

    when:

    script:
    """
    samtools view -bS $sam -o ${sample}_markers.bam
    rm $sam
    """
}
