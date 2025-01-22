process prefetch_and_split {
    tag "$sra_id"

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_*.fastq.gz")

    publishDir "${params.readsdir}", mode: 'copy'

    script:
    """
    echo "Prefetching $sra_id to S3"
    prefetch $sra_id -v -v

    echo "Splitting $sra_id into FASTQ files directly on S3"
    fasterq-dump $sra_id --split-files --threads ${task.cpus} -v -v

    echo "Compressing FASTQ files"
    gzip ${sra_id}_1.fastq
    gzip ${sra_id}_2.fastq
    """
}