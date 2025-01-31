process FILTER_AND_GLUE_CONTIGS {

    tag "$sample_id"
    label "process_low"
    // container 'pegi3s/biopython:1.78'
    container 'ilepeli/adar:0.0.2'

    input: 
    val(sample_id)
    path(mummer_tile_file)
    path(chosen_ref)
    path(ref_ch)
    path(mummer_delta_file)

    output:
    path("*intermediate_scaffold.fasta"), emit: scaffold


    script:

    """

    reference_fasta=\$(cat ${chosen_ref})
    ref_name=\${reference_fasta%.*}

    filter_and_glue_contigs.py \\
    ${mummer_tile_file} \\
    \$reference_fasta \\
    ${mummer_delta_file} \\
    --out_scaffold_name "${sample_id}_\${ref_name}_intermediate_scaffold.fasta"

    """

}
