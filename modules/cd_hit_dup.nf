process CD_HIT_DUP {

    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/cd-hit-auxtools:4.8.1--h4ac6f70_3'
    maxRetries 0

    input: 
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_dedup.fastq.gz'), emit: reads

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // note: using zcat because cd-hit-dup can't intake gzipped reads

    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        cd-hit-dup \\
            -i ${prefix}.fastq.gz \\
            -o ${prefix}_dedup.fastq.gz \\
        """
    } else {

        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

        ## cd-hit-dup seg faults if using with -i, -i2, so just call it for r1 and r2

        cd-hit-dup \\
            -i <(zcat ${prefix}_1.fastq.gz) -o ${prefix}_1_dedup.fastq.gz

        cd-hit-dup \\
            -i <(zcat ${prefix}_2.fastq.gz) -o ${prefix}_2_dedup.fastq.gz
        """
    }

}
