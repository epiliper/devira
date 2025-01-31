process CHOOSE_BEST_REF {

    tag "$sample_id"
    label 'process_high'
    container 'staphb/skani:0.2.2'

    input: 
    tuple val(sample_id), path(contigs)
    path(ref_fastas)

    output: 
    tuple val(sample_id), path(contigs), path("./CHOSEN_REF"), emit: chosen_ref

    script:

    // note: misc skani args taken from 
    // https://github.com/broadinstitute/viral-assemble/blob/master/assembly.py
    """
    skani dist ${contigs} ${ref_fastas} \\
    -m 50 \\
    -s 50 \\
    -c 20 \\
    --min-af 15 \\
    --no-learned-ani \\
    --robust \\
    --detailed \\
    --ci \\
    -n 10 \\
    --no-marker-index \\
    -o ${sample_id}.refs_skani_dist.full.tsv

    CHOSEN_REF_FASTA=\$(cut -f 1 "${sample_id}.refs_skani_dist.full.tsv" | tail +2 | head -1)

    echo "\$CHOSEN_REF_FASTA" > CHOSEN_REF

    """
}
