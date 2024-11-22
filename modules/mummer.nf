process MUMMER {
    tag "$sample_id"
    label 'process_high'
    container 'ilepeli/adar:0.0.2'

    input: 
    val(sample_id)
    path(contigs_fasta)
    path(chosen_ref)
    path(ref_ch)

    output:
    path("${sample_id}_post_filter.delta"), emit: delta
    path("${sample_id}.tiling"), emit: tile

    script:

    """
    chosen_ref=\$(cat ${chosen_ref})

    touch ${sample_id}_post_filter.delta

    nucmer --prefix ${sample_id}_pre_filter --maxgap 200 --minmatch 10 \$chosen_ref ${contigs_fasta}

    delta-filter ${sample_id}_pre_filter.delta > ${sample_id}_post_filter.delta

    show-tiling -a -i 60.0 -l 200 -v 30.0 ${sample_id}_post_filter.delta > ${sample_id}.tiling

    """


    
}
