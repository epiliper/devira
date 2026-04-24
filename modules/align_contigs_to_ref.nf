process ALIGN_CONTIGS_TO_REF {
    tag "${task.ext.prefix}"
    label 'process_single'
    container "quay.io/biocontainers/bwa:0.7.19--h577a1d6_1"

    input: 
    tuple val(meta), val(tax_info), path(contigs), val(ref_info), path(ref)

    output:

    tuple val(meta), val(tax_info), val(ref_info), path("*.sam"), emit: alignment

    script:
    def prefix = task.ext.prefix
    """
    echo ${prefix}.sam
    bwa index $ref
    bwa mem -B 1 -k 15 -a $ref $contigs > ${prefix}.sam
    """
}
