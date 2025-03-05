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

    ch_reads
        .cross(MAKE_REFERENCE_FASTA.out.ref)
        .multiMap { it -> 
            reads: it[0]
            ref: it[1]
        }
        .set { ch_output }

    emit:
    contigs = contigs
    reads   = ch_output.reads
    ref     = ch_output.ref
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
