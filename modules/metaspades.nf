process METASPADES {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/spades:4.0.0'


    input:
    tuple val(meta), path(reads)


    output:
    path("contigs.fasta"), emit: contigs_fasta
    path("corrected/*.fastq.gz"), emit: reads_corrected
    val("${meta.id}"), emit: sample_id

    script: 

    """

    if [ "${meta.single_end}" = true ]; then
        metaspades.py \\
        -k 21 \\
        -s ${reads} \\
        -o .

    else
        metaspades.py \\
        -k 21 \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o .
    fi
    """

}
