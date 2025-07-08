process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata", mode: 'link'
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    //errorStrategy 'retry'
    //maxRetries 3

    input:
    tuple val(sample), path(reads)

    output:
    val(sample), emit: sample
    path("${sample}_kneaddata.fastq.gz"), emit: fastq
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    shell:
    
    """
    echo $sample

    kneaddata --unpaired $reads \
              --reference-db ${params.human_genome} --output ./ \
              --processes ${task.cpus} --output-prefix ${sample}_kneaddata \
              --trimmomatic /cluster/tufts/bonhamlab/shared/conda-envs/metaphlan_v4.2/.CondaPkg/.pixi/envs/default/share/trimmomatic

    gzip ${sample}_kneaddata*.fastq
    """  
}
