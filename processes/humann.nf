

process humann {
    // process samples with humann
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
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
    path "${sample}_0.log"
    path "${sample}_3_reactions.tsv"
    path "${sample}_4_pathabundance.tsv"

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \
        --protein-database ${params.humann_db}/humann_protein_db \
        --nucleotide-database ${params.humann_nucleotide_db} \
        --utility-database ${params.humann_utility_db} \
        --output-basename $sample 
    """
}

process humann_regroup {
    // regroup gene families into functional categories
    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"

    input:
    val  sample
    path genefamilies

    output:
    val  sample , emit: sample
    path "${sample}_ecs.tsv", emit: ecs
    path "${sample}_kos.tsv", emit: kos
    path "${sample}_pfams.tsv", emit: pfams


    script:

    """
    humann_regroup_table --input $genefamilies --output ${sample}_ecs.tsv --custom ${params.humann_utility_db}/map_level4ec_uniclust90.txt.gz
    humann_regroup_table --input $genefamilies --output ${sample}_kos.tsv --custom ${params.humann_utility_db}/map_ko_uniref90.txt.gz
    humann_regroup_table --input $genefamilies --output ${sample}_pfams.tsv --custom ${params.humann_utility_db}/map_pfam_uniref90.txt.gz
    """
}   


process humann_rename {
    // add names to individual feature IDs
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"

    input:
    val sample
    path ecs
    path kos
    path pfams

    output:
    val  sample , emit: sample
    path "${sample}_ecs_rename.tsv"
    path "${sample}_kos_rename.tsv"
    path "${sample}_pfams_rename.tsv"

    script:

    // Rename file output to include humann DB used for functional profiling

    """
    humann_rename_table --input $ecs   --output ${sample}_ecs_${params.humann_version}_rename.tsv   --custom ${params.humann_utility_db}/map_level4ec_name.txt.gz
    humann_rename_table --input $kos   --output ${sample}_kos_${params.humann_version}_rename.tsv   --custom ${params.humann_utility_db}/map_ko_name.txt.gz
    humann_rename_table --input $pfams --output ${sample}_pfams_${params.humann_version}_rename.tsv --custom ${params.humann_utility_db}/map_pfam_name.txt.gz
    """
}

