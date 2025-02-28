process GAPFILL_WITH_REF {
    tag "$meta.id"
    label 'process_low'
    container 'ilepeli/adar:0.0.2'

    input:
    tuple val(meta), path(intermediate_contigs), path(chosen_ref), val(ref_name), path(reads)

    output:
    tuple val(meta), path(reads), path("*_imputed.fasta"), val(ref_name), emit: prep_scaffold

    script:

    def prefix = "${meta.id}_${ref_name}"
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
