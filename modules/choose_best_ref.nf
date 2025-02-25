process CHOOSE_BEST_REF {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/skani:0.2.2'

    input: 
    tuple val(meta), path(contigs)
    path(ref_fastas)

    output: 
    tuple val(meta), path(contigs), path("*_ref.fasta"), emit: chosen_ref


    script:
    def prefix = task.ext.prefix ?: ''

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
    -o ${prefix}.refs_skani_dist.full.tsv

    CHOSEN_REF_FASTA=\$(cut -f 1 "${prefix}.refs_skani_dist.full.tsv" | tail +2 | head -1)
    cp \$CHOSEN_REF_FASTA \${CHOSEN_REF_FASTA%*.}_ref.fasta

    """
}
