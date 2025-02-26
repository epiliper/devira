process GAPFILL_WITH_REF {
    tag "$meta.id"
    label 'process_low'
    container 'ilepeli/adar:0.0.2'

    input:
    tuple val(meta), path(intermediate_contigs), path(chosen_ref)

    output:
    path("*_imputed.fasta"), emit: prep_scaffolds

    script:
    """
    ref_name=\$(cat $chosen_ref | grep \\> | tr -d \\> | cut -d ' ' -f 1)
    prefix="${meta.id}_\${ref_name}"

    nucmer $chosen_ref ${intermediate_contigs} --prefix \${prefix}
    delta-filter \${prefix}.delta > \${prefix}_filtered.delta

    prep_scaffolds.py \\
        \${prefix}_filtered.delta \\
        $intermediate_contigs \\
        $chosen_ref \\
        \${prefix}_imputed.fasta
    """
}
