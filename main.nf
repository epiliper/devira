//  check args
if (params.input) { ch_input = file(params.input) } else { exit 1, 'No samplesheet specified!' }
if (!params.refs) { exit 1, "Reference genome multifasta not specified!" }
if (params.run_kraken2 & params.kraken2_db == null ) {exit 1, "Must provide a kraken2 database with --kraken2_db!" }

// Import subworkflows
include { INPUT_CHECK               } from './subworkflows/input_check'
include { FASTP_MULTIQC             } from './subworkflows/fastp_multiqc'
include { CONTIG_GEN                } from './subworkflows/contig_gen'
include { PROFILE_READS             } from './subworkflows/profile_reads'

// Import modules
include { KRAKEN2                   } from './modules/kraken2'
include { CD_HIT_DUP                } from './modules/cd_hit_dup'
include { SUBSAMPLE_FASTQ           } from './modules/subsample'
include { CHOOSE_BEST_REF           } from './modules/choose_best_ref'
include { ORDER_AND_ORIENT          } from './modules/order_and_orient'
include { MUMMER                    } from './modules/mummer'
include { FILTER_AND_GLUE_CONTIGS   } from './modules/filter_and_glue_contigs'
include { GAPFILL_GAP2SEQ           } from './modules/gapfill_gap2seq'

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

    PROFILE_READS(
        ch_sample_input,
        file(params.kraken2_db),
        file(params.taxids)
    )

    ch_sample_input = PROFILE_READS.out.profiled_reads

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

    MUMMER(
        CHOOSE_BEST_REF.out.chosen_ref
    )

    FILTER_AND_GLUE_CONTIGS(
        MUMMER.out.delta_tile
    )

    GAPFILL_GAP2SEQ(
        FILTER_AND_GLUE_CONTIGS
        .out
        .intermediate_scaffold
        .join(ch_sample_input, by: [0])
    )
}
