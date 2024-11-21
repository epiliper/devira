process IMPUTE_FROM_REFERENCE {

    tag "${sample_id}"
    label 'process_high'
    container 'ilepeli/viral_assemble:latest'

    input: 
    val(sample_id)
    path(chosen_ref)
    path(gapfilled_fasta)
    path(refs)

    output:
    path("${sample_id}.scaffolded_imputed.fasta"), emit: scaffolded_imputed_fasta

    script:
    
    """
    chosen_ref_fasta=\$(cat ${chosen_ref})

    assembly.py impute_from_reference \\
    ${gapfilled_fasta} \\
    \$chosen_ref_fasta \\
    ${sample_id}.scaffolded_imputed.fasta \\
    --newName ${sample_id} \\
    --loglevel=DEBUG

    """
}
