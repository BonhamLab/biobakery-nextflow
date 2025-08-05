process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", mode: 'link', pattern: "*.tsv"
    publishDir "$params.outdir/metaphlan", mode: 'link', pattern: "*.sam" 
    
    input:
    val(sample)
    path(kneads)
    path(paired), optional true
    path(unpaired), optional true

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_bowtie2.tsv" , emit: bowtie2
    path "${sample}.sam"         , emit: sam
    path "${sample}_concatenated.fastq.gz", emit: concatenated

    script:
    // metaphlan4 changed metaphlan db variable from bowtie2db to db_dir
    // also changed from bowtie2out to mapout

    if (params.paired_end == false){
        kneads = preprocess_paired_end_metaphlan(sample, paired, unpaired)
    }

    mp_ver = params.metaphlan_version
    if (mp_ver == 'metaphlan_v4') {
        db_arg = 'db_dir'
        out_arg = 'mapout';
    }else if (mp_ver == 'metaphlan_v3'){
        db_arg = 'bowtie2db'
        out_arg = 'bowtie2out'
    }else {
        throw new Exception("The metaphlan_version must be 'metaphlan_v4' or 'metaphlan_v3', got '${mp_ver}'")
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
 