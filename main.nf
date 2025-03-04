//  check args
if (params.input) { ch_input = file(params.input) } else { exit 1, 'No samplesheet specified!' }
if (!params.refs) { exit 1, "Reference genome multifasta not specified!" }
if (params.run_kraken2 & params.kraken2_db == null ) {exit 1, "Must provide a kraken2 database with --kraken2_db!" }

// Import subworkflows
include { INPUT_CHECK       } from './subworkflows/input_check'
include { FASTP_MULTIQC     } from './subworkflows/fastp_multiqc'
include { CONTIG_GEN        } from './subworkflows/contig_gen'
include { PROFILE_READS     } from './subworkflows/profile_reads'

// Import modules
include { KRAKEN2           } from './modules/kraken2'
include { CD_HIT_DUP        } from './modules/cd_hit_dup'
include { SUBSAMPLE_FASTQ   } from './modules/subsample'
include { BWA_MEM2_ALIGN    } from './modules/bwa_align'
include { HAPLOTYPE_CALLER  } from './modules/haplotype_caller'
include { CALL_CONSENSUS    } from './modules/call_consensus'

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

    //ref_ch = Channel.fromPath("${params.refs}/*.fasta").collect()

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
        file(params.refs)
    )

    BWA_MEM2_ALIGN(
        CONTIG_GEN.out.contigs,
    )

    HAPLOTYPE_CALLER(
        BWA_MEM2_ALIGN.out.bam
        .join(BWA_MEM2_ALIGN.out.ref)
    )

    CALL_CONSENSUS(
        HAPLOTYPE_CALLER.out.vcf
        .join(HAPLOTYPE_CALLER.out.gatk_deps)
    )
}
