//  check args
if (params.input) { ch_input = file(params.input) } else { exit 1, 'No samplesheet specified!' }
if (!params.refs || !file(params.refs).exists()) { exit 1, "Invalid path specified: reference genome multifasta (--ref) " }
if (!params.kraken2_db || !file(params.kraken2_db).exists()) {exit 1, "Invalid path specified: kraken2 database (--kraken2_db)" }
if (!params.taxids || !file(params.taxids).exists()) {exit 1, "Invalid path specified: taxon ID TSV (--taxids)" }

// Import subworkflows
include { INPUT_CHECK           } from './subworkflows/input_check'
include { FASTP_MULTIQC         } from './subworkflows/fastp_multiqc'
include { PROFILE_READS         } from './subworkflows/profile_reads'
include { CONTIG_GEN            } from './subworkflows/contig_gen'
include { CONSENSUS_ASSEMBLY    } from './subworkflows/consensus_assembly'

// Import modules
include { KRAKEN2           } from './modules/kraken2'
include { SUBSAMPLE_FASTQ   } from './modules/subsample'
include { SUMMARY           } from './modules/summary'
                                                              
log.info("`7MM===Yb.        `7MMF'   `7MF'db  `7MM===Mq.        db      ")
log.info("  MM    `Yb.        `MA     ,V        MM   `MM.      ;MM:     ")
log.info("  MM     `Mb  .gP=Ya VM:   ,V `7MM    MM   ,M9      ,V^MM.    ")
log.info("  MM      MM ,M'   Yb MM.  M'   MM    MMmmdM9      ,M  `MM    ")
log.info("  MM     ,MP 8M====== `MM A'    MM    MM  YM.      AbmmmqMA   ")
log.info("  MM    ,dP' YM.    ,  :MM;     MM    MM   `Mb.   A'     VML  ")
log.info(".JMMmmmdP'    `Mbmmd'   VF    .JMML..JMML. .JMM..AMA.   .AMMA.\n\n")
                                                              
workflow {

    INPUT_CHECK (
        ch_input
    )

    inreads = INPUT_CHECK.out.reads

    if (params.subsample_raw) {

        SUBSAMPLE_FASTQ(
            inreads,
            params.subsample_raw
        )

        inreads = SUBSAMPLE_FASTQ.out.reads
    }

    FASTP_MULTIQC (
        inreads,
        params.adapter_fasta,
        params.save_merged,
        params.skip_fastp,
        params.skip_multiqc
    )

    ch_sample_input = FASTP_MULTIQC.out.reads

    PROFILE_READS(
        ch_sample_input,
        file(params.kraken2_db),
        file(params.taxids)
    )

    ch_sample_input = PROFILE_READS.out.profiled_reads

    CONTIG_GEN(
        ch_sample_input,
        params.contig_method,
        file(params.refs),
        params.min_contig_length
    )

    CONSENSUS_ASSEMBLY(
        CONTIG_GEN.out.contigs,
        CONTIG_GEN.out.reads,
    )

    // [ meta, tax_info, log, contigs, n50 ]
    FASTP_MULTIQC.out.trim_log
    .combine(
        //CONTIG_GEN.out.contigs.map { meta, tax_info, ref_info, scaffold -> 
        //[ meta, tax_info, scaffold ] }
        //.join(CONTIG_GEN.out.contig_stats, by: [0, 1]),
        //by: 0)
        CONTIG_GEN.out.contigs,
        by: 0)
    .map { meta, log, tax_info, ref_info, contigs  -> [ meta, tax_info, ref_info, log, contigs ] }
    .set { contig_sum_ch }

    // [ meta, tax_info, log, contigs, n50, consensus, init_covstats, final_covstats]
    ch_summary_in = contig_sum_ch
    .join(
        CONTIG_GEN.out.contig_stats
            .join(CONSENSUS_ASSEMBLY.out.final_consensus, by: [0, 1, 2])
            .join(CONSENSUS_ASSEMBLY.out.init_covstats, by: [0, 1, 2])
            .join(CONSENSUS_ASSEMBLY.out.final_covstats, by: [0, 1, 2]),
    by: [0, 1, 2])
    .map { meta, tax_info, ref_info, log, contigs, contig_stats, consensus, init_covstats, final_covstats -> 
    [ meta, tax_info, ref_info, log, contigs, contig_stats, consensus, init_covstats, final_covstats ] }

    SUMMARY(
        ch_summary_in
    )

    SUMMARY.out.summary
    .collectFile(storeDir: "${params.output}", name:"${params.run_name}_summary.tsv", keepHeader: true, sort: {file -> file.name })

}
