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
    def temp_out = "${meta.id}_temp.fasta"
    def consensus = "${prefix}_consensus.fasta"


    """
    gatk IndexFeatureFile -I $vcf
    gatk FastaAlternateReferenceMaker \\
        -R $ref \\
        -O $temp_out \\
        -V $vcf

    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' $temp_out > $consensus
    """
}
