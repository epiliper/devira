include { SKANI                 } from '../modules/skani'
include { SELECT_REFERENCES     } from '../modules/select_references'
include { MAKE_REFERENCE_FASTA  } from '../modules/make_reference_fasta'

workflow REFERENCE_PREP {
    take: 
    contigs  // path(meta, contigs)
    ch_reads // path(meta, reads)
    db       // path(db)

    main:
    SKANI(contigs, db)

    SELECT_REFERENCES(SKANI.out.dist)

    SELECT_REFERENCES.out.refs_tsv
    .map { meta, refs_tsv -> add_ref_info_to_meta(meta, refs_tsv) }
    .flatten().collate( 2,false )
    .set { ch_new_meta }

    MAKE_REFERENCE_FASTA(
        ch_new_meta,
        db
    )

    // Combine reads with references while keeping meta
    ch_reads
        .combine(MAKE_REFERENCE_FASTA.out.ref, by: [0])
        .set { ch_reads_refs }

    // Combine with contigs (since contigs also have meta)
    ch_reads_refs
        .combine(contigs, by: [0])
        .multiMap { it ->
            reads: [ it[0], it[1] ]
            ref: [ it[0], it[2], it[3] ]
            contigs: [ it[0], it[4] ]
        }
        .set { ch_output }

    emit:
    reads       = ch_output.reads
    ref         = ch_output.ref 
    contigs     = ch_output.contigs
}

def add_ref_info_to_meta(meta, refs_tsv) {
    def new_meta_list = []

    refs_tsv.splitCsv(sep: '\t', skip: 1, header: false).collect { row ->
        def new_meta = []
        new_meta = [ meta, [ acc: row[0], tag: row[1], header: row[2].toString() ] ]
        new_meta_list.add(new_meta)
    }
    
    return new_meta_list
}
