process FILL_GAPS_WITH_REFERENCE {
    tag "$sample_id"
    label 'process_low'
    container 'ilepeli/adar:0.0.2'

    input:
    tuple (
            val(sample_id),
            path(intermediate_contigs),
            path(chosen_ref)
          )

    output:
    path("*imputed_assembly.fasta"), emit: prep_scaffolds

    script:

    """
    reference_fasta=$chosen_ref
    ref_name=\${reference_fasta%%.*}
    prefix=${sample_id}_\$ref_name

    nucmer $chosen_ref $intermediate_contigs --prefix \${prefix}
    delta-filter \${prefix}.delta > \${prefix}_filtered.delta

    prep_scaffolds.py \\
        \${prefix}_filtered.delta \\
        $intermediate_contigs \\
        \$reference_fasta \\
        \${prefix}_imputed_assembly.fasta

    """

}
