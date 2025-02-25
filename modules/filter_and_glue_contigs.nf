process FILTER_AND_GLUE_CONTIGS {

    tag "$meta.id"
    label "process_low"
    container 'quay.io/epil02/adar:0.0.4'

    input: 

    tuple val(meta), path(mummer_delta_file), path(mummer_tile_file), path(chosen_ref)

    output:
    tuple val(meta), path("*intermediate_scaffold.fasta"), path(chosen_ref), emit: intermediate_scaffold

    script:

    def prefix = "${meta.id}"

    """

    reference_fasta=${chosen_ref}
    ref_name=\${reference_fasta%.*}

    filter_and_glue_contigs.py \\
    ${mummer_tile_file} \\
    $chosen_ref \\
    ${mummer_delta_file} \\
    --out_scaffold_name "${prefix}_\${ref_name}_intermediate_scaffold.fasta"

    """

}
