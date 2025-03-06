process MAKE_REFERENCE_FASTA {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"                                                              
    label 'process_single'                                                         
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    tuple val(meta), val(ref_info)
    path db

    output:
    tuple val(meta), val(ref_info), path("*.fa"), emit: ref

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    samtools faidx ${db} ${ref_info.acc} > ${prefix}.fa
    """
}

