//  check args
if (params.input) { ch_input = file(params.input) } else { exit 1, 'No samplesheet specified!' }
// if (!params.refs) { exit 1, "Reference genome multifasta not specified!" }
if (params.run_kraken2 & params.kraken2_db == null ) {exit 1, "Must provide a kraken2 database with --kraken2_db!" }

include { INPUT_CHECK } from './subworkflows/input_check'
include { KRAKEN2 } from './modules/kraken2'
include { FASTP_MULTIQC } from './subworkflows/fastp_multiqc'
include { CD_HIT_DUP } from './modules/cd_hit_dup'


log.info "Reference-guided denovo assembly"

workflow {
    log.info("Workflow started")

    INPUT_CHECK (
        ch_input
    )

    FASTP_MULTIQC (
        INPUT_CHECK.out.reads,
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
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
        )

        if (params.kraken2_variants_host_filter) {
            ch_sample_input = KRAKEN2.out.unclassified_reads_fastq
        }
    }

    if (!params.skip_dedup) {

        CD_HIT_DUP(
            ch_sample_input
        )


    }
}
