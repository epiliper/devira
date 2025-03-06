process SELECT_REFERENCES {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/epil02/revica-strm:0.0.4'

    input: 
    tuple val(meta), path(contigs), path(dist_report)

    output:
    tuple val(meta), path("*failed_assembly.tsv"), optional: true, emit: failed_assembly
    tuple val(meta), path("*covstats.tsv"), optional: true, emit: refs_tsv

    script: 
    """
    select_references.py --sample ${meta.id} \\
        --dist_report $dist_report \\
        --ani_thres 60.0 \\
        --align_ref_thres 60.0
    """
}
