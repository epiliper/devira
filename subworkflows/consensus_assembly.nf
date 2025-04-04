include { BWA_MEM2_ALIGN as BWA_MEM2_INITIAL } from '../modules/bwa_align'
include { BWA_MEM2_ALIGN as BWA_MEM2_FINAL   } from '../modules/bwa_align'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_INITIAL } from '../modules/ivar_consensus'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_FINAL   } from '../modules/ivar_consensus'

workflow CONSENSUS_ASSEMBLY {
    take:
    contigs
    reads

    main:
    BWA_MEM2_INITIAL(
        contigs
        .join(reads, by: [0, 1])
    )

    IVAR_CONSENSUS_INITIAL(
        BWA_MEM2_INITIAL.out.bam
        .join(BWA_MEM2_INITIAL.out.ref, by: [0, 1])
    )

    BWA_MEM2_FINAL(
        IVAR_CONSENSUS_INITIAL.out.consensus
        .join(reads, by: [0, 1])
    )

    IVAR_CONSENSUS_FINAL(
        BWA_MEM2_FINAL.out.bam
        .join(BWA_MEM2_FINAL.out.ref, by: [0, 1])
    )

    emit:
    init_alignments = BWA_MEM2_INITIAL.out.bam
    init_consensus = IVAR_CONSENSUS_INITIAL.out.consensus
    init_covstats = BWA_MEM2_INITIAL.out.covstats

    final_alignments = BWA_MEM2_FINAL.out.bam
    final_consensus = IVAR_CONSENSUS_FINAL.out.consensus
    final_covstats = BWA_MEM2_FINAL.out.covstats
}
