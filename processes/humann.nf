

process humann {
    // process samples with humann
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main", mode: 'link'
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }

    errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    val  sample
    path catkneads
    path profile

    output:
    val  sample                       , emit: sample
    path "${sample}_2_genefamilies.tsv" , emit: genefamilies
    path "${sample}_0.log",             emit: log
    path "${sample}_3_reactions.tsv",   emit: reactions
    path "${sample}_4_pathabundance.tsv", emit: pathabundance

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \
        --protein-database ${params.humann_db}/chocophlan \
        --nucleotide-database ${params.humann_nucleotide_db}/uniref \
        --utility-database ${params.humann_utility_db}/utility_mapping \
        --output-basename $sample 
    """
}

process humann_rename {
    // rename humann output 
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"

    input:
    val sample
    path genefamilies
    path log
    path reactions
    path pathabundance

    output:
    val  sample , emit: sample
    path "${sample}_2_genefamilies.tsv"
    path "${sample}_0.log"    
    path "${sample}_3_reactions.tsv"
    path "${sample}_4_pathabundance.tsv"

    script:
    hp_ver = params.humann_version
    if (hp_ver == 'humann_v37') {
        huamnn_db_ver = 'v201901' // need to double check this value
    }else if (hp_ver == 'humann_v4a'){
        huamnn_db_ver = 'v201901_v31'

    }else {
        throw new Exception("The humann_version must be 'humann_v37' or 'humann_v4a', got '${hp_ver}'")
    }

    // Rename file output to include humann DB used for functional profiling

    """
    mv "${sample}_2_genefamilies.tsv" "${sample}_2_genefamilies_${huamnn_db_ver}.tsv"
    mv "${sample}_0.log"  "${sample}_0_${huamnn_db_ver}.log"
    mv "${sample}_3_reactions.tsv"  "${sample}_3_reactions_${huamnn_db_ver}.tsv"
    mv "${sample}_4_pathabundance.tsv"  "${sample}_4_pathabundanc_${huamnn_db_ver}.tsv"
    """
}

