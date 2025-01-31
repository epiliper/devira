process FILTER_AND_GLUE_CONTIGS {

    tag "$sample_id"
    label "process_low"
    // container 'pegi3s/biopython:1.78'
    container 'ilepeli/adar:0.0.2'

    input: 

    tuple val(sample_id), path(mummer_delta_file), path(mummer_tile_file), path(chosen_ref)

    output:
    tuple val(sample_id), path("*intermediate_scaffold.fasta"), path(chosen_ref), emit: intermediate_scaffold


    script:

    """

    reference_fasta=${chosen_ref}
    ref_name=\${reference_fasta%.*}

    filter_and_glue_contigs.py \\
    ${mummer_tile_file} \\
    $chosen_ref \\
    ${mummer_delta_file} \\
    --out_scaffold_name "${sample_id}_\${ref_name}_intermediate_scaffold.fasta"

    """

}
