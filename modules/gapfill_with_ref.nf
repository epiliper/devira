process GAPFILL_WITH_REF {
    tag "${task.ext.prefix}"
    label 'process_low'
    container 'ilepeli/adar:0.0.2'

    input:
    tuple val(meta), val(tax_info), path(intermediate_contigs), path(chosen_ref), val(ref_info), path(reads)

    output:
    tuple val(meta), val(tax_info), val(ref_info), path("*_imputed.fasta"), emit: prep_scaffold

    script:

    def prefix = task.ext.prefix
    """
    nucmer $chosen_ref ${intermediate_contigs} --prefix ${prefix}
    delta-filter ${prefix}.delta > ${prefix}_filtered.delta

    prep_scaffolds.py \\
        ${prefix}_filtered.delta \\
        $intermediate_contigs \\
        $chosen_ref \\
        ${prefix}_imputed.fasta
    """
}
