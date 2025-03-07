process MUMMER {
    tag "${task.ext.prefix}"
    label 'process_high'
    container 'quay.io/epil02/adar:0.0.4'

    input: 
    tuple val(meta), val(tax_info), val(ref_info), path(ref), path(contigs)

    output:

    tuple val(meta), val(tax_info), path("*_post_filter.delta"), path("*.tiling"), path(ref), val(ref_info), emit: delta_tile

    script:
    def prefix = task.ext.prefix

    """
    sample_id=${meta[0]}

    touch ${prefix}_post_filter.delta

    nucmer --prefix ${prefix}_pre_filter --maxgap 200 --minmatch 10 $ref ${contigs}

    delta-filter ${prefix}_pre_filter.delta > ${prefix}_post_filter.delta

    show-tiling -a -i 60.0 -l 200 -v 30.0 ${prefix}_post_filter.delta > ${prefix}.tiling

    """
}
