process bam2fastq {
    tag "bam2fastq $sample"

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fastq")

    shell:
    
    """
    echo $sample

    samtools fastq -@ ${task.cpus} ${reads} > ${sample}.fastq
    """  
}
