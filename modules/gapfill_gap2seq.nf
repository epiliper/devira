process GAPFILL_GAP2SEQ {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.4'

    input: 
    tuple val(meta), path(scaffold_fasta), path(chosen_ref), path(reads)

    output:
    tuple val(meta), path("*gapfilled.fasta"), path(reads), emit: gapfilled_scaffold

    script:

    """
    scaffold_fasta_path=${scaffold_fasta}
    out_scaffold=\${scaffold_fasta_path%*.}_gapfilled.fasta

    gapfill.py \\
        --in_scaffold $scaffold_fasta \\
        --in_reads $reads \\
        --out_scaffold \$out_scaffold
    """
}
