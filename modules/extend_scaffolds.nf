process EXTEND_SCAFFOLDS {
    tag "${task.ext.prefix}"
    label "process_high"
    container "quay.io/jefffurlong/biopython_biocode:latest"

    input:
    tuple val(meta), val(tax_info), path(intermediate_scaffold), path(chosen_ref), val(ref_info), path(contigs)

    output:
    tuple val(meta), val(tax_info), path("*scaffold_extended.fa"), path(chosen_ref), val(ref_info), emit: extended_scaffold

    script:

    def out = "${task.ext.prefix}_scaffold_extended.fa"
    def stats_out ="${task.ext.prefix}_contig_stats.tsv"

    """
    extend_gaps.py \\
        -s $intermediate_scaffold \\
        -r $chosen_ref \\
        -q $contigs \\
        -o $out \\
        -t $task.cpus
    """
}
