process CALL_CONSENSUS {
    tag "${meta.id}_${ref_name}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(ref), val(ref_name), path(vcf), path(fai), path(dict)

    output:
    tuple val(meta), path("*_consensus.fasta"), emit: consensus
    
    script:
    def prefix = "${meta.id}_${ref_name}"
    def consensus = "${prefix}_consensus.fasta"


    """
    gatk IndexFeatureFile -I $vcf
    gatk FastaAlternateReferenceMaker \\
        -R $ref \\
        -O $consensus \\
        -V $vcf
    """
}
