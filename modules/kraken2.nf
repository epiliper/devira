process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0'

    // this workflow has been modified to classify contigs, not reads.
    input:
    tuple val(meta), path(contigs)
    path db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_contigs
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_contigs
    tuple val(meta), path('*classifiedcontigs.txt')   , optional:true, emit: classified_contigs_assignment
    tuple val(meta), path('*report.txt')                           , emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def classified   = "${prefix}.classified.fasta"
    def unclassified = "${prefix}.unclassified.fasta"
    def classified_option           = "--classified-out ${classified}"
    def unclassified_option         = "--unclassified-out ${unclassified}"
    def readclassification_option   = "--output ${prefix}.kraken2.classifiedcontigs.txt"

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $args \\
        $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
