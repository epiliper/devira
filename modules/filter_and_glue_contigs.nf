process FILTER_AND_GLUE_CONTIGS {
    tag "${task.ext.prefix}"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.4'

    input: 
    tuple val(meta), val(tax_info), val(ref_info), path(ref), path(contigs)

    output:

    tuple val(meta), val(tax_info), val(ref_info), path("*intermediate_scaffold.fasta"), emit: intermediate_scaffold
    tuple val(meta), val(tax_info), val(ref_info), path("*_contig_stats.tsv"), emit: contig_stats

    script:
    def prefix = task.ext.prefix
    def tiling = "${prefix}.tiling"
    def post_filter = "${prefix}_post_filter.delta"
    def pre_filter = "${prefix}_pre_filter"
    def contig_stats = "${prefix}_contig_stats.tsv"

    """
    touch $post_filter
    nucmer --prefix $pre_filter --maxgap 200 --minmatch 10 $ref $contigs
    delta-filter ${pre_filter}.delta > $post_filter
    show-tiling -a -i 60.0 -l 200 -v 30.0 $post_filter > $tiling

    filter_and_glue_contigs.py \\
        $tiling \\
        $ref \\
        $post_filter \\
        --out_scaffold_name "${prefix}_intermediate_scaffold.fasta"
    
    ## get contig stats
    contig_stats.py -c $contigs \\
        -t $tiling \\
        -o $contig_stats
    """
}
