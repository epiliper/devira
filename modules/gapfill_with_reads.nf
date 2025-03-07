process GAPFILL_WITH_READS {
    tag "${task.ext.prefix}"
    label 'process_medium'
    container 'quay.io/epil02/adar:0.0.6'

    input: 
    tuple val(meta), val(tax_info), path(scaffold_fasta), path(chosen_ref), val(ref_info), path(reads)

    output:
    tuple val(meta), val(tax_info), path("*gapfilled.fasta"), path(chosen_ref), val(ref_info), path(reads), emit: gapfilled_scaffold

    script:
    def reads_in = meta.single_end ? "${reads}" : "${reads[0]},${reads[1]}"
    def prefix = task.ext.prefix

    """
    scaffold_fasta_path=${scaffold_fasta}
    out_scaffold=${prefix}_gapfilled_temp.fasta
    final_scaffold=${prefix}_gapfilled.fasta

    gapfill.py \\
        --in_scaffold $scaffold_fasta \\
        --in_reads $reads_in \\
        --out_scaffold \$out_scaffold

    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' \$out_scaffold > \$final_scaffold
    rm \$out_scaffold
    """
}
