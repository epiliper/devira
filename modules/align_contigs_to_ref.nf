process ALIGN_CONTIGS_TO_REF {
    tag "${task.ext.prefix}"
    label 'process_single'
    container "quay.io/epil02/revica-strm:0.0.4"

    input: 
    tuple val(meta), val(tax_info), path(contigs), val(ref_info), path(ref)

    output:

    tuple val(meta), val(tax_info), val(ref_info), path("*.sam"), emit: alignment

    script:
    def prefix = task.ext.prefix
    """
    echo ${prefix}.sam
    bwa-mem2 index $ref
    bwa-mem2 mem -B 1 -k 15-a $ref $contigs > ${prefix}.sam
    """
}
