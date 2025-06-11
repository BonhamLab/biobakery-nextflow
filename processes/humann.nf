

process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }

    errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    val  sample
    path profile
    path catkneads

    output:
    val  sample                       , emit: sample
    path "${sample}_genefamilies.tsv" , emit: genefamilies
    path "${sample}_pathabundance.tsv"
    path "${sample}_pathcoverage.tsv"

    script:

    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output \ # add --search-mode uniref90 for 3.7
        --output-basename $sample \
        --protein-database ${params.humann_protein_db} \
        --nucleotide-database ${params.humann_nucleotide_db} \
        --utility-mapping ${params.humann_utility_db} \
        --metaphlan-options="--index mpa_vOct22_CHOCOPhlAnSGB_202403 --bowtie2db ${params.metaphlan_db} -t rel_ab_with_read_stats"

    """
}

process humann_regroup {
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
    humann_config --update database_folders utility_mapping `realpath ${params.humann_utility_db}`
    humann_regroup_table --input $genefamilies --output ${sample}_ecs.tsv --groups uniref90_level4ec
    humann_regroup_table --input $genefamilies --output ${sample}_kos.tsv --groups uniref90_ko
    humann_regroup_table --input $genefamilies --output ${sample}_pfams.tsv --groups uniref90_pfam
    """
}   

process humann_rename {
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

    """
    humann_config --update database_folders utility_mapping `realpath ${params.humann_utility_db}`
    humann_rename_table --input $ecs   --output ${sample}_ecs_rename.tsv   --names ec
    humann_rename_table --input $kos   --output ${sample}_kos_rename.tsv   --names kegg-orthology
    humann_rename_table --input $pfams --output ${sample}_pfams_rename.tsv --names pfam
    """
}