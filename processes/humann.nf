

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
    path "${sample}_genefamilies_${params.humann_version}.tsv" , emit: genefamilies
    path "${sample}_${params.humann_version}.log", optional:true           
    path "${sample}_reactions_${params.humann_version}.tsv"
    path "${sample}_pathabundance_${params.humann_version}.tsv"
    path "${sample}_pathcoverage_${params.humann_version}.tsv", optional:true

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \
        --protein-database ${params.humann_db}/uniref \
        --nucleotide-database ${params.humann_db}/chocophlan \
        --utility-database ${params.humann_db}/utility_mapping \
        --output-basename $sample 

    if [[ "$params.humann_version" == 'humann_v4a' ]]; then 
        mv "${sample}_2_genefamilies.tsv" "${sample}_genefamilies_${params.humann_version}.tsv" 
        mv "${sample}_0.log" "${sample}_${params.humann_version}.log"
        mv "${sample}_3_reactions.tsv" "${sample}_reactions_${params.humann_version}.tsv"
        mv "${sample}_4_pathabundance.tsv" "${sample}_pathabundance_${params.humann_version}.tsv"

    elif [[ "$params.humann_version" == 'humann_v37' ]]; then 

        mv "${sample}_genefamilies.tsv" "${sample}_genefamilies_${params.humann_version}.tsv" 
        mv "${sample}_pathabundance.tsv" "${sample}_pathabundance_${params.humann_version}.tsv"
        mv "${sample}_pathcoverage.tsv" "${sample}_pathcoverage_${params.humann_version}.tsv"
    
    fi
    """
}
