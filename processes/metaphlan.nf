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
    // metaphlan4 changed metaphlan db variable from bowtie2db to db_dir
    // also changed from bowtie2out to mapout
    if (params.metaphlan_version == 'metaphlan_v4') {
    db_arg = 'db_dir'
    out_arg = 'mapout'}
    else (params.metaphlan_version == 'metaphlan_v3'){
    db_arg = 'bowtie2db'
    out_arg = 'bowtie2out'
    }
    
    // Run metaphlan on samples
    """
    metaphlan $kneads -o ${sample}_profile.tsv \
        --${out_arg} ${sample}_bowtie2.tsv \
        --samout ${sample}.sam \
        --input_type fastq \
        --nproc ${task.cpus} \
        --${db_arg} ${params.metaphlan_db} \
        --index ${params.metaphlan_index} \
        -t rel_ab_w_read_stats
    """
}

    process rename_metaphlan_database_version {
    // Rename file output to include metaphlan DB used for taxonomic ID
    input:
    val(sample)

    output:
    path "${sample}_profile_${params.metaphlan_index}.tsv"
    path "${sample}_bowtie2_${params.metaphlan_index}.tsv"
    path "${sample}_${params.metaphlan_index}.sam"

    script:
    """
    mv "${sample}_profile.tsv" "${sample}_profile_${params.metaphlan_index}.tsv"
    mv "${sample}_bowtie2.tsv" "${sample}_bowtie2_${params.metaphlan_index}.tsv"
    mv "${sample}.sam"  "${sample}_${params.metaphlan_index}.sam"
    """
    }
 