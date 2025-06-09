process bam2fastq {
    tag "bam2fastq $sample"
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("{sample}.fastq")

    shell:
    
    """
    echo $sample

    samtools fastq -@ {task.cpus} > {sample}.fastq
    """  
}
