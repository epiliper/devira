process CONTIG_STATS {
    tag "${task.ext.prefix}"
    label "process_low"
    container "quay.io/jefffurlong/biopython_biocode:latest"

    input: 
    tuple val(meta), val(tax_info), path(mummer_delta), path(mummer_tile), path(chosen_ref), val(ref_info), path(contigs)

    output:
    tuple val(meta), val(tax_info), val(ref_info), path("*_contig_stats.tsv"), emit: contig_stats

    script:
    def prefix = task.ext.prefix
    def out = "${prefix}_contig_stats.tsv"

    """
    contig_stats.py -c $contigs \\
        -t $mummer_tile \\
        -o $out
    """
}
