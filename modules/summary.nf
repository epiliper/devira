process SUMMARY {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_single'
    container 'quay.io/epil02/revica-strm:0.0.4'

    input:
    tuple val(meta), val(tax_info), val(ref_info), path(fastp_trim_log), path(contigs), path(contig_stats), path(consensus), path(init_covstats), path(final_covstats)

    output:
    path("*.tsv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''

    """
    taxon=${tax_info.name}
    taxon_reads=${tax_info.num_reads}
    contig_header=\$(sed -n '1p' $contig_stats | tr '\t\n\r' '   ')
    contig_values=\$(sed -n '2p' $contig_stats | tr '\t\n\r' '   ')

    # raw reads and trimmed reads
    raw_reads=\$(grep -A1 "before filtering:" ${fastp_trim_log} | grep 'total reads:' | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}')
    trimmed_reads=\$(grep -A1 "after filtering:" ${fastp_trim_log} | grep 'total reads:' | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}')
    pct_reads_trimmed=\$(echo "\${trimmed_reads}/\${raw_reads}*100" | bc -l)
    pct_reads_trimmed_formatted=\$(printf "%.2f" "\${pct_reads_trimmed}")

    #####################
    # ALIGN TO SCAFFOLD #
    #####################

    # mapped reads
    #mapped_reads_ref=\$(awk -F'\t' '{print \$NF}' $init_covstats)
    mapped_reads_ref=\$(awk 'NR==2 {print \$NF}' $init_covstats | tr -d ' \r\n')
    pct_reads_mapped=\$(echo "\${mapped_reads_ref}/\${raw_reads}*100" | bc -l)
    pct_reads_mapped_formatted_ref=\$(printf "%.2f" "\${pct_reads_mapped}")

    # coverage/depth
    coverage_ref=\$(awk 'BEGIN {FS="\t"} NR==2 {print \$5}' $init_covstats)
    mean_depth_ref=\$(awk 'BEGIN {FS="\t"} NR==2 {print \$6}' $init_covstats)

    ##############################
    # ALIGN TO INITIAL CONSENSUS #
    ##############################

    # mapped reads
    #mapped_reads_c1=\$(awk -F'\t' '{print \$NF}' $final_covstats)
    mapped_reads_c1=\$(awk 'NR==2 {print \$NF}' $final_covstats | tr -d ' \r\n')
    pct_reads_mapped=\$(echo "\${mapped_reads_c1}/\${raw_reads}*100" | bc -l)
    pct_reads_mapped_formatted_c1=\$(printf "%.2f" "\${pct_reads_mapped}")

    # coverage/depth
    coverage_c1=\$(awk 'BEGIN {FS="\t"} NR==2 {print \$5}' $final_covstats)
    mean_depth_c1=\$(awk 'BEGIN {FS="\t"} NR==2 {print \$6}' $final_covstats)

    ####################
    # CONSENSUS GENOME #
    ####################

    consensus_length=\$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${consensus} | awk 'FNR==2{print val,\$1}')
    num_ns_consensus=\$(grep -v "^>" ${consensus} | tr -c -d N | wc -c)
    pct_ns=\$(echo "\${num_ns_consensus}/\${consensus_length}*100" | bc -l | awk 'FNR==1{print val,\$1}')
    pct_ns_formatted=\$(printf "%.4f" "\${pct_ns}")
    num_as_consensus=\$(grep -v "^>" ${consensus} | tr -c -d A | wc -c)
    num_cs_consensus=\$(grep -v "^>" ${consensus} | tr -c -d C | wc -c)
    num_gs_consensus=\$(grep -v "^>" ${consensus} | tr -c -d G | wc -c)
    num_ts_consensus=\$(grep -v "^>" ${consensus} | tr -c -d T | wc -c)
    num_non_ns_ambiguous=\$(echo "\${consensus_length}-\${num_as_consensus}-\${num_cs_consensus}-\${num_gs_consensus}-\${num_ts_consensus}-\${num_ns_consensus}" | bc -l)
    
    ##################
    # OUTPUT TO TSV #
    ##################

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "sample_name" "raw_reads" "trimmed_reads" "pct_reads_trimmed" "taxon" "taxon_reads"\
        \$contig_header "ref_tag" "ref_acc"\
        "mapped_reads_ref" "pct_reads_mapped_ref" "coverage_ref" "mean_depth_ref" \
        "mapped_reads_c1" "pct_reads_mapped_c1" "coverage_c1" "mean_depth_c1" \
        "consensus_length" "num_ns_consensus" "pct_ns" "num_ambiguous" > ${prefix}_summary.tsv

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$meta.id" "\$raw_reads" "\$trimmed_reads" "\$pct_reads_trimmed_formatted" "\$taxon" "\$taxon_reads" \
        \$contig_values "${ref_info.tag}" "${ref_info.acc}"\
        "\$mapped_reads_ref" "\$pct_reads_mapped_formatted_ref" "\$coverage_ref" "\$mean_depth_ref" \
        "\$mapped_reads_c1" "\$pct_reads_mapped_formatted_c1" "\$coverage_c1" "\$mean_depth_c1" \
        "\$consensus_length" "\$num_ns_consensus" "\$pct_ns_formatted" "\$num_non_ns_ambiguous" \
        >> "${prefix}_summary.tsv"
    """
}
