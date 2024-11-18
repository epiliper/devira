process METASPADES {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/spades:4.0.0'


    input:
    tuple val(meta), path(reads)


    output:
    path("contigs.fasta"), emit: contigs_fasta
    path("corrected/*.fastq.gz"), emit: reads_corrected

    script: 

    """

    if [ "${meta.single_end}" = true ]; then
        metaspades.py \\
        -s ${reads} \\
        -o .
    else
        metaspades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o .
    fi
    """

}
