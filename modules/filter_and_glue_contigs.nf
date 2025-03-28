process FILTER_AND_GLUE_CONTIGS {
    tag "${task.ext.prefix}"
    label 'process_single'
    container "quay.io/epil02/gap2seq_bio_pysam:0.0.2"


    input: 
    tuple val(meta), val(tax_info), val(ref_info), path(aln), path(contigs), path(reads)

    output:
    tuple val(meta), val(tax_info), val(ref_info), path("*intermediate_scaffold.fasta"), emit: intermediate_scaffold
    tuple val(meta), val(tax_info), val(ref_info), path("*_contig_stats.tsv"), emit: contig_stats

    script:
    def prefix = task.ext.prefix
    def contig_stats = "${prefix}_contig_stats.tsv"
    def reads_in = meta.single_end ? "$reads" : "${reads[0]},${reads[1]}"

    """
    scaffold.py \\
        -a $aln \\
        -q $contigs \\
        -r $reads_in \\
        -o ${prefix}_intermediate_scaffold.fasta \\
        -p ${prefix} -c $contig_stats
    """
}
