process HAPLOTYPE_CALLER {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_name)

    output:
    tuple val(meta), path("*.gcvf"), emit: gcvf

    script:
    def prefix = "${meta.id}_${ref_name}"
    def bam_rg = "${prefix}_rg.bam"
    """
    samtools faidx $ref
    gatk CreateSequenceDictionary -R $ref

    gatk HaplotypeCaller -I $bam \\
        -R $ref -O ${prefix}_calls.gcvf -A AlleleFraction \\
        --native-pair-hmm-threads ${task.cpus - 1} \\
    """
}
