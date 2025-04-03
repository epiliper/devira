process SKANI {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/skani:0.2.2'

    input: 
    tuple val(meta), val(tax_info), path(contigs)
    path(ref)

    output: 
    tuple val(meta), val(tax_info), path(contigs), path("*dist.full.tsv"), emit: dist

    script:
    def prefix = "$meta.id"

    // note: misc skani args taken from 
    // https://github.com/broadinstitute/viral-assemble/blob/master/assembly.py
    """
    skani dist ${contigs} ${ref} \\
    -m 50 \\
    -s 50 \\
    -c 20 \\
    --min-af 60 \\
    --no-learned-ani \\
    --robust \\
    --detailed \\
    --ci \\
    --ri \\
    --no-marker-index \\
    -o ${prefix}.refs_skani_dist.full.tsv
    """
}
