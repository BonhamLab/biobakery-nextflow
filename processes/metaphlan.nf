process metaphlan {
    tag "metaphlan on $sample"
    // publishDir "$params.outdir/metaphlan", pattern: "*.tsv" // once fix for sam compression is found
    publishDir "$params.outdir/metaphlan" // keeps sam file

    input:
    val(sample)
    path(kneads)

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_bowtie2.tsv" , emit: bowtie2
    path "${sample}.sam"         , emit: sam


    script:
    // metphlan4 changed metaphlan db variable from bowtie2db to db_dir
    // also changed from bowtie2out to mapout
    if (params.metaphaln_ver == 'metaphlan4') {
    db_arg = 'db_dir'
    out_arg = 'mapout'}
    else (params.metaphaln_ver == 'metaphlan3.1.0'){
    db_arg = 'bowtie2db'
    out_arg = 'bowtie2out'
    }
    
    """
    metaphlan $kneads -o ${sample}_profile.tsv \
        --${out_arg} ${sample}_bowtie2.tsv \
        --samout ${sample}.sam \
        --input_type fastq \
        --nproc ${task.cpus} \
        --${dbarg} ${params.metaphlan_db} \
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
    # Duplicate headers in output - skipping header validation
    samtools view -bS --no-PG ${sam} -o ${sample}_markers.bam
 
    # Alternative approach: strip and rebuild header
    # samtools view -S ${sam} | samtools view -b -o ${sample}_markers.bam
    
    rm $sam
    """
}
