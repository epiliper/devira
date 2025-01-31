//  check args
if (params.input) { ch_input = file(params.input) } else { exit 1, 'No samplesheet specified!' }
if (!params.refs) { exit 1, "Reference genome multifasta not specified!" }
if (params.run_kraken2 & params.kraken2_db == null ) {exit 1, "Must provide a kraken2 database with --kraken2_db!" }

include { INPUT_CHECK } from './subworkflows/input_check'
include { KRAKEN2 } from './modules/kraken2'
include { FASTP_MULTIQC } from './subworkflows/fastp_multiqc'
include { CD_HIT_DUP } from './modules/cd_hit_dup'
include { SUBSAMPLE_FASTQ } from './modules/subsample'

include { CONTIG_GEN } from './subworkflows/contig_gen'

include { CHOOSE_BEST_REF } from './modules/choose_best_ref'
include { ORDER_AND_ORIENT } from './modules/order_and_orient'
include { GAP_FILL } from './modules/gap_fill'
include { IMPUTE_FROM_REFERENCE } from './modules/impute_from_reference'
include { MUMMER } from './modules/mummer'
include { FILTER_AND_GLUE_CONTIGS } from './modules/filter_and_glue_contigs'
include { PREP_SCAFFOLDS } from './modules/prep_scaffolds'



log.info("   █████████   ██████████     █████████   ███████████      ") 
log.info("  ███░░░░░███ ░░███░░░░███   ███░░░░░███ ░░███░░░░░███     ") 
log.info(" ░███    ░███  ░███   ░░███ ░███    ░███  ░███    ░███     ") 
log.info(" ░███████████  ░███    ░███ ░███████████  ░██████████      ") 
log.info(" ░███░░░░░███  ░███    ░███ ░███░░░░░███  ░███░░░░░███     ")
log.info(" ░███    ░███  ░███    ███  ░███    ░███  ░███    ░███     ")
log.info(" █████   █████ ██████████   █████   █████ █████   █████    ")
log.info("░░░░░   ░░░░░ ░░░░░░░░░░   ░░░░░   ░░░░░ ░░░░░   ░░░░░     ") 
log.info "Assisted de-novo assembly via reference\n"


workflow {

    INPUT_CHECK (
        ch_input
    )

    ref_ch = Channel.fromPath("${params.refs}/*.fasta").collect()

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

    if (params.run_kraken2) {
        KRAKEN2 (
            ch_sample_input,
            file(params.kraken2_db),
            true,
            false
        )

        ch_sample_input = KRAKEN2.out.classified_reads_fastq
    }

    if (!params.skip_dedup) {

        CD_HIT_DUP(
            ch_sample_input
        )
        ch_sample_input = CD_HIT_DUP.out.reads
    }

    CONTIG_GEN(
        ch_sample_input,
        params.contig_method,
    )

    CHOOSE_BEST_REF(
        CONTIG_GEN.out.contigs,
        ref_ch
    )

    // ORDER_AND_ORIENT(
    MUMMER(
        CHOOSE_BEST_REF.out.chosen_ref,
        ref_ch
    )

    FILTER_AND_GLUE_CONTIGS(
        MUMMER.out.delta_tile,
        ref_ch
        // CONTIG_GEN.out.sid,
        // MUMMER.out.tile,
        // CHOOSE_BEST_REF.out.chosen_ref,
        // ref_ch,
        // MUMMER.out.delta,
    )

    PREP_SCAFFOLDS(
        FILTER_AND_GLUE_CONTIGS
        .out
        .intermediate_scaffold,
        ref_ch
        // CONTIG_GEN.out.sid,
        // FILTER_AND_GLUE_CONTIGS.out.scaffold,
        // CHOOSE_BEST_REF.out.chosen_ref,
    )



    // GAP_FILL(
    //     CONTIG_GEN.out.sid,
    //     ch_sample_input,
    //     ORDER_AND_ORIENT.out.intermediate_scaffold
    // )

    // IMPUTE_FROM_REFERENCE(
    //     CONTIG_GEN.out.sid,
    //     CHOOSE_BEST_REF.out.chosen_ref,
    //     GAP_FILL.out.gapfilled_fasta,
    //     ref_ch,
    // )
}
