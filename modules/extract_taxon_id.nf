process EXTRACT_TAXON_ID {
    tag "$meta.id - $taxid"
    label 'process_single'
    container 'quay.io/epil02/adar:0.0.4'

    input:
    tuple val(meta), path(classified_fastq), path(kraken_output), path(kraken2_report), val(taxid)

    output:
    tuple val(meta), path("*_TAX.fastq.gz"), emit: tax_reads
    path("*_profile.tsv"), emit: profile_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_in = meta.single_end ? "-s1 ${classified_fastq}" : "-s1 ${classified_fastq[0]} -s2 ${classified_fastq[1]}"
    def reads_out = meta.single_end ? "-o ${prefix}_${taxid}_TAX.fastq" : "-o ${prefix}_${taxid}_TAX.fastq -o2 ${prefix}_${taxid}_2_TAX.fastq"

    """
    extract_kraken_reads.py \\
        -k $kraken_output \\
        -r $kraken2_report \\
        $reads_in \\
        -t $taxid \\
        $reads_out \\
        --fastq-output \\
        --include-children \\
        --include-parents > /dev/null

    gzip *.fastq

    numlines=\$(zcat ${prefix}_${taxid}_TAX.fastq.gz | wc -l)
    numreads=\$(echo "\$numlines / 4" | bc)

    if [ "\$(echo "\$numreads >= 2000" | bc)" -eq 0 ]; then
        echo "Insufficient reads (\$numreads), skipping..."
        rm *.fastq.gz
        touch FAILED_TAX.fastq.gz
    else
        echo "Read threshold met: found \$numreads reads"
    fi

    printf "sample\ttaxon_id\tnum_found_reads\n" > ${meta.id}_${taxid}_profile.tsv
    printf "${meta.id}\t$taxid\t\$numreads\n" >> ${meta.id}_${taxid}_profile.tsv
    """
}
