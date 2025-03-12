include { KRAKEN2 } from '../modules/kraken2'
include { EXTRACT_TAXON_ID } from '../modules/extract_taxon_id'

// this workflow classifies fastq files (ideally QCed and trimmed beforehand)
// against a kraken2 database, after which the classified reads are grouped into files 
// associated with a specific taxonID

// extract_kraken_reads is naively called on each taxonID specified in .tsv file
// with output fastqs only making it through the pipeline if they match or exceed 
// a given threshold read count.

workflow PROFILE_CONTIGS {

    take:
    contigs_meta  // tuple val(meta), path(contigs)
    krak_db     // path to Kraken2 database
    taxids      // path to taxonid list

    main:
    KRAKEN2(contigs_meta, krak_db, true, false)

    // combine into a tuple for convenience
    KRAKEN2.out.classified_contigs
        .join(KRAKEN2.out.classified_contigs_assignment)
        .join(KRAKEN2.out.report)
        .set { krak_out_ch } 

    // for each taxonID, duplicate the classified reads and attempt to 
    // extract matching reads
    Channel.fromPath(taxids)
        .splitCsv(header: false, strip: true, sep: '\t')
        .map { row -> row[0] }
        .set { tid_ch }

    krak_out_ch
        .combine(tid_ch)
        .set { extract_input_ch }

    EXTRACT_TAXON_ID(extract_input_ch)

    // filter taxonID fastqs to have a significant number of reads
    EXTRACT_TAXON_ID
        .out.tax_contigs
        .filter { meta, taxid, num_bases_f, contigs -> contigs.name != "FAILED_TAX.fasta" }
        .map { meta, taxid, num_bases_f, contigs -> [ meta, [ taxid: taxid, num_bases: num_bases_f.text.toInteger() ], contigs ]}
        .set { profiled_contigs }

    EXTRACT_TAXON_ID.out.profile_report
    .collectFile(
        storeDir: "${params.output}", 
        name: "${params.run_name}_profile.tsv",
        keepHeader: true,
        sort: {file -> file.text}
    )

    emit:
    profiled_contigs          = profiled_contigs
}
