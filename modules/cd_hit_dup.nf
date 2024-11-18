process CD_HIT_DUP {

    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/cd-hit-auxtools:4.8.1--h4ac6f70_3'
    maxRetries 0

    input: 
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_dedup.fastq'), emit: reads

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

        gunzip ${prefix}_1.fastq.gz 
        gunzip ${prefix}_2.fastq.gz

        ## cd-hit-dup seg faults if using with -i, -i2 with gzipped input, so unzip first

        cd-hit-dup \\
            -i ${prefix}_1.fastq -o ${prefix}_1_dedup.fastq \\
            -i2 ${prefix}_2.fastq -o2 ${prefix}_2_dedup.fastq

        """
    }

}
