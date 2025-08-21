

process humann {
    // process samples with humann
    tag "humann on $sample"
    publishDir "$params.outdir/humann/$params.humann_version", mode: 'link'
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }

    errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    val(sample)
    path(catkneads)
    path(profile)

    output:
    val  sample                       , emit: sample
    path "${sample}_2_genefamilies.tsv" , emit: genefamilies
    path "${sample}_0.log",             emit: log
    path "${sample}_3_reactions.tsv",   emit: reactions
    path "${sample}_4_pathabundance.tsv", optional:true, emit: pathabundance

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \
        --protein-database ${params.humann_db}/uniref \
        --nucleotide-database ${params.humann_db}/chocophlan \
        --utility-database ${params.humann_db}/utility_mapping \
        --output-basename $sample 
    """
}

process humann_rename {
    // rename humann output 
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/$params.humann_version"

    input:
    val(sample)
    path(genefamilies)
    path(log)       
    path(reactions)  
    path(pathabundance), optional true

    output:
    path "${sample}_2_genefamilies_${params.humann_version}.tsv"
    path "${sample}_0_${params.humann_version}.log"
    path "${sample}_3_reactions_${params.humann_version}.tsv"
    path "${sample}_4_pathabundance_${params.humann_version}.tsv", optional:true
    
    script:
    // Rename file output to include humann DB used for functional profiling
    """
    hp_ver="${params.humann_version}"
    mv $genefamilies "${sample}_2_genefamilies_${hp_ver}.tsv"
    mv $log  "${sample}_0_${hp_ver}.log"
    mv $reactions  "${sample}_3_reactions_${hp_ver}.tsv"

    if [[ "$hp_ver" == "humann_v4a" ]]; then
        mv $pathabundance  "${sample}_4_pathabundance_${hp_ver}.tsv"
    fi
    """
}
