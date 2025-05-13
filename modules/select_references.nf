process SELECT_REFERENCES {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/epil02/revica-strm:0.0.4'

    input: 
    tuple val(meta), val(tax_info), path(contigs), path(dist_report)
    path(refs)

    output:
    tuple val(meta), val(tax_info), path("*failed_assembly.tsv"), optional: true, emit: failed_assembly
    tuple val(meta), val(tax_info), path("*covstats.tsv"), optional: true, emit: refs_tsv

    script: 
    """
    select_references.py --sample ${meta.id} \\
        --dist_report $dist_report \\
        --ani_thres 60.0 \\
        --align_ref_thres 50.0 \\
        --reference_fasta $refs
    """
}
