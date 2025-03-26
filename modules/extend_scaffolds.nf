process EXTEND_SCAFFOLDS {
    tag "${task.ext.prefix}"
    label "process_high"
    container "quay.io/jefffurlong/biopython_biocode:latest"

    input:
    tuple val(meta), val(tax_info), val(ref_info), path(intermediate_scaffold), path(query_fasta)

    output:
    tuple val(meta), val(tax_info), val(ref_info), path("*scaffold_extended.fa"), emit: extended_scaffold

    script:

    def out = "${task.ext.prefix}_scaffold_extended.fa"
    def stats_out ="${task.ext.prefix}_contig_stats.tsv"
    def args = task.ext.args ?: ''

    """
    extend_gaps.py \\
        -s $intermediate_scaffold \\
        -q $query_fasta \\
        -o $out \\
        -t $task.cpus \\
        $args
    """
}
