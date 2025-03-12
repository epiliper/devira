process EXTRACT_TAXON_ID {
    tag "$meta.id - $taxid"
    label 'process_single'
    container 'quay.io/epil02/adar:0.0.4'

    input:
    tuple val(meta), path(classified_contigs), path(kraken_output), path(kraken2_report), val(taxid)

    output:
    tuple val(meta), val(taxid), path("NUM_BASES"), path("*_TAX.fasta"), emit: tax_contigs
    path("*_profile.tsv"), emit: profile_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def contigs_out = "${prefix}_${taxid}_TAX.fasta"

    """
    extract_kraken_reads.py \\
        -k $kraken_output \\
        -r $kraken2_report \\
        -s $classified_contigs \\
        -t $taxid \\
        -o $contigs_out \\
        --include-children \\
        --include-parents > /dev/null

    ## count number of characters that aren't part of the header
    num_bases=\$(cat $contigs_out | grep -v \\> | wc -m | tr -d ' ')

    if [ "\$(echo "\$num_bases >= 500" | bc)" -eq 0 ]; then
        echo "Insufficient bases (\$num_bases), skipping..."
        rm *.fasta
        touch FAILED_TAX.fasta
    else
        echo "Base count threshold met: found \$num_bases bases"
    fi

    echo \$num_bases > NUM_BASES

    printf "sample\ttaxon_id\tnum_associated_bases\n" > ${meta.id}_${taxid}_profile.tsv
    printf "${meta.id}\t$taxid\t\$num_bases\n" >> ${meta.id}_${taxid}_profile.tsv
    """
}
