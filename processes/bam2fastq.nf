process bam2fastq {
    tag "bam2fastq $sample"
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)

    output:

    shell:
    
    """
    echo $sample

    kneaddata --input ${reads[0]} --input ${reads[1]} \
              --reference-db $human_genome --output ./ \
              --processes ${task.cpus} --output-prefix ${sample}_kneaddata \
              --trimmomatic /opt/conda/share/trimmomatic

    gzip *.fastq
    """  
}
