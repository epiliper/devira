process GAP_FILL {

    tag "${sample_id}"
    label "process_high"
    container 'ilepeli/viral_assemble:latest'


    input: 
    val(sample_id)
    tuple val(meta), path(reads) 
    path(intermediate_scaffold)

    output:
    path("./${sample_id}.intermediate_gapfill.fasta"), emit: gapfilled_fasta

    script:

    def read_files_comma_sep = reads.join(', ')
    log.info("${read_files_comma_sep}")

    if (meta.single_end) {
        """
        samtools import -@ ${task.cpus} -0 ${reads[0]} -o ${sample_id}.bam
        """
    } else {
        """
        samtools import -@ ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}.bam
        """
    }

    """

    assembly.py gapfill_gap2seq \\
    ${intermediate_scaffold} \\
    ${sample_id}.bam \\
    ${sample_id}.intermediate_gapfill.fasta \\
    --memLimitGb ${task.mem} \\
    --maskErrors \\
    --loglevel=DEBUG

    """
}
