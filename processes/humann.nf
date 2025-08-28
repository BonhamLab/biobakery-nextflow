

process humann {
    // process samples with humann
    tag "humann on $sample"
    publishDir "$params.outdir/humann/$params.humann_version", mode: 'link'
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }

    // errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
    // maxRetries 3


    input:
    val(sample)
    path(catkneads)
    path(profile)

    output:
    val  sample                       , emit: sample
    path "${sample}_2_genefamilies_${params.humann_version}.tsv" , emit: genefamilies
    path "${sample}_0_${params.humann_version}.log",             emit: humann_log
    path "${sample}_3_reactions_${params.humann_version}.tsv",   emit: reactions
    path "${sample}_4_pathabundance_${params.humann_version}.tsv", emit: pathabundance

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \
        --protein-database ${params.humann_db}/uniref \
        --nucleotide-database ${params.humann_db}/chocophlan \
        --utility-database ${params.humann_db}/utility_mapping \
        --output-basename $sample 

    mv "${sample}_2_genefamilies.tsv" "${sample}_2_genefamilies_${params.humann_version}.tsv" 
    mv "${sample}_0.log" "${sample}_0_${params.humann_version}.log"
    mv "${sample}_3_reactions.tsv" "${sample}_3_reactions_${params.humann_version}.tsv"
    mv "${sample}_4_pathabundance.tsv" "${sample}_4_pathabundance_${params.humann_version}.tsv"
    """
}
