process MUMMER {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.4'

    input: 
    tuple val(meta), path(contigs_fasta), path(chosen_ref)

    output:

    tuple val(meta), path("*_post_filter.delta"), path("*.tiling"), path(chosen_ref), emit: delta_tile

    script:
    def prefix = "${meta.id}"

    """
    sample_id=${meta[0]}

    touch ${prefix}_post_filter.delta

    nucmer --prefix ${prefix}_pre_filter --maxgap 200 --minmatch 10 $chosen_ref ${contigs_fasta}

    delta-filter ${prefix}_pre_filter.delta > ${prefix}_post_filter.delta

    show-tiling -a -i 60.0 -l 200 -v 30.0 ${prefix}_post_filter.delta > ${prefix}.tiling

    """


    
}
