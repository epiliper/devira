process CHOOSE_BEST_REF {

    tag "$fasta"
    label 'process_high'
    container 'ilepeli/viral_assemble:latest'

    input: 
    val(sample_id)
    path(contigs)
    path(ref_fastas)

    output: 
    path("./CHOSEN_REF_BASENAME"), emit: chosen_ref

    script:
    """

    assembly.py skani_contigs_to_refs \\
    ${contigs}  \\
    ${ref_fastas} \\
    ${sample_id}.refs_skani_dist.full.tsv \\
    ${sample_id}.refs_skani_dist.top.tsv \\
    ${sample_id}.ref_clusters.tsv \\
    --loglevel=DEBUG


    CHOSEN_REF_FASTA=\$(cut -f 1 "${sample_id}.refs_skani_dist.full.tsv" | tail +2 | head -1)
          
    basename "\$CHOSEN_REF_FASTA" .fasta > CHOSEN_REF_BASENAME

    """
}
