process GAPFILL_WITH_READS {
    tag "${meta.id}_${ref_name}"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.5'

    input: 
    tuple val(meta), path(scaffold_fasta), path(chosen_ref), val(ref_name), path(reads)

    output:
    tuple val(meta), path("*gapfilled.fasta"), path(chosen_ref), val(ref_name), path(reads), emit: gapfilled_scaffold

    script:
    def fastq_to_bam = meta.single_end ? "-0 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def prefix = "${meta.id}_${ref_name}"
    def bam = "${prefix}.bam"

    """
    samtools import $fastq_to_bam -o $bam
    scaffold_fasta_path=${scaffold_fasta}
    out_scaffold=${prefix}_gapfilled_temp.fasta
    final_scaffold=${prefix}_gapfilled.fasta

    gapfill.py \\
        --in_scaffold $scaffold_fasta \\
        --in_reads $bam \\
        --out_scaffold \$out_scaffold

    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' \$out_scaffold > \$final_scaffold
    rm \$out_scaffold
    """
}
