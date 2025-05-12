include { SKANI                 } from '../modules/skani'
include { SELECT_REFERENCES     } from '../modules/select_references'
include { MAKE_REFERENCE_FASTA  } from '../modules/make_reference_fasta'

workflow REFERENCE_PREP {
    take: 
    contigs  // path(meta, tax_info, contigs)
    ch_reads // path(meta, tax_info, reads)
    db       // path(db)

    main:
    SKANI(contigs, db)

    SELECT_REFERENCES(SKANI.out.dist, db)

    // if we have the same reference being selected for multiple taxon ids, then select
    // only the taxon id with the highest read count to be used to generate said reference
    // this should not be a problem with kraken2 databases that are low redundancy
    SELECT_REFERENCES.out.refs_tsv
    .map { meta, tax_info, refs_tsv -> add_ref_info_to_meta(meta, tax_info, refs_tsv) }
    .flatten().collate( 3,false )
    .collect(flat: false, sort: { a ,b -> b[1].num_reads <=> a[1].num_reads } )
    .flatten().collate( 3,false)
    .unique { meta, tax_info, ref_info -> [ meta.id, ref_info.tag ]}
    .set { ch_new_meta }

    MAKE_REFERENCE_FASTA(
        ch_new_meta,
        db
    )

    // Combine reads with references while keeping meta
    // [ meta, tax_info, ref_info, ref, reads ]
    MAKE_REFERENCE_FASTA.out.ref
        .combine(ch_reads, by: [0, 1])
        .set { ch_reads_refs }

    // Combine with contigs (since contigs also have meta)
    // [ meta, tax_info, ref_info, ref, reads, contigs ]
    ch_reads_refs
        .combine(contigs, by: [0, 1])
        .multiMap { it ->
            reads:      [ it[0], it[1], it[4] ]
            ref:        [ it[0], it[1], it[2], it[3] ]
            contigs:    [ it[0], it[1], it[5] ]
        }
        .set { ch_output }

    emit:
    reads       = ch_output.reads
    ref         = ch_output.ref 
    contigs     = ch_output.contigs
}

def add_ref_info_to_meta(meta, tax_info, refs_tsv) {
    def new_meta_list = []

    refs_tsv.splitCsv(sep: '\t', skip: 1, header: false).collect { row ->
        def new_meta = []
        new_meta = [ meta, tax_info, [ acc: row[0], tag: row[1], header: row[2].toString() ] ]
        new_meta_list.add(new_meta)
    }
    
    return new_meta_list
}
