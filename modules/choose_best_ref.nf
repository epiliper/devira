process CHOOSE_BEST_REF {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/skani:0.2.2'

    input: 
    tuple val(meta), path(contigs)
    path(ref)

    output: 
    tuple val(meta), path(contigs), path("*_ref.fasta"),                 emit: chosen_ref
    tuple val(meta), path("REF_ACC"), path("REF_DESC"), path("REF_TAG"), emit: ref_info

    script:
    def prefix = task.ext.prefix ?: ''

    // note: misc skani args taken from 
    // https://github.com/broadinstitute/viral-assemble/blob/master/assembly.py
    """

    skani dist ${contigs} ${ref} \\
    -m 50 \\
    -s 50 \\
    -c 20 \\
    --min-af 15 \\
    --no-learned-ani \\
    --robust \\
    --detailed \\
    --ci \\
    -n 10 \\
    --ri \\
    --no-marker-index \\
    -o ${prefix}.refs_skani_dist.full.tsv

    # gather ref info and store in file
    ref_acc=\$(cut -f 6 "${prefix}.refs_skani_dist.full.tsv" | tail +2 | head -1| cut -d ' ' -f 1)
    ref_tag=\$(cut -f 6 "${prefix}.refs_skani_dist.full.tsv" | tail +2 | head -1| cut -d ' ' -f 2)
    ref_desc=\$(cut -f 6 "${prefix}.refs_skani_dist.full.tsv" | tail +2 | head -1| cut -d ' ' -f 3-)
    echo \$ref_tag > REF_TAG
    echo \$ref_desc > REF_DESC
    echo \$ref_acc > REF_ACC

    # make single fasta from selected ref for downstream work
    cat $ref | grep \$ref_acc -A 1 > "\${ref_tag}_ref.fasta"

    """
}
