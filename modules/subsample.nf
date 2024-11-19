process SUBSAMPLE_FASTQ {

    tag "$meta.id"
    label 'process_high'
    container 'staphb/rasusa:2.1.0'

    input: 
    tuple val(meta), path(reads)
    val (n_reads)

    output: 
    tuple val(meta), path("*_sub.fastq.gz"), emit: reads

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def SEED = "100"

    if (meta.single_end) {

        def out = "${prefix}_sub.fastq.gz"
    
        """

        rasusa reads -s ${SEED} \\
        --num ${n_reads} \\
        -o ${out} \\
        ${reads} 

        """
    } else {

        def out1 = "${prefix}_1_sub.fastq.gz"
        def out2 = "${prefix}_2_sub.fastq.gz"

        """

        rasusa reads -s ${SEED} \\
        --num ${n_reads} \\
        -o ${out1} \\
        -o ${out2} \\
        ${reads[0]} \\
        ${reads[1]}

        """
    }
}
