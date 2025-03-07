process MEGAHIT {

    tag "$meta.id"
    label 'process_high'
    container 'vout/megahit:release-v1.2.9'


    input:
    tuple val(meta), val(tax_info), path(reads)


    output:
    tuple val(meta), val(tax_info), path("${meta.id}/final.contigs.fa"), emit: contigs

    script: 

    """

    if [ "${meta.single_end}" = true ]; then
        megahit \\
        -r ${reads} \\
        --presets meta-sensitive \\
        -o ${meta.id}

    else
        megahit \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${meta.id}
    fi
    """

}
