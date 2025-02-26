process GAPFILL_WITH_READS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.5'

    input: 
    tuple val(meta), path(scaffold_fasta), path(chosen_ref), path(reads)

    output:
    tuple val(meta), path("*gapfilled.fasta"), path(chosen_ref), path(reads), emit: gapfilled_scaffold

    script:
    def fastq_to_bam = meta.single_end ? "-0 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def bam = "${meta.id}.bam"

    """
    samtools import $fastq_to_bam -o $bam
    scaffold_fasta_path=${scaffold_fasta}
    out_scaffold=\${scaffold_fasta_path%*.}_gapfilled.fasta

    gapfill.py \\
        --in_scaffold $scaffold_fasta \\
        --in_reads $bam \\
        --out_scaffold \$out_scaffold
    """
}
