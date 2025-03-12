include { MUMMER                    } from '../modules/mummer'
include { FILTER_AND_GLUE_CONTIGS   } from '../modules/filter_and_glue_contigs'
include { GAPFILL_WITH_READS        } from '../modules/gapfill_with_reads'
include { GAPFILL_WITH_REF          } from '../modules/gapfill_with_ref'
include { REFERENCE_PREP            } from './reference_prep'

workflow REFINE_CONTIGS {

   take: 
   contigs_meta
   reads_meta       // tuple val(meta), path(reads)
   ref_ch           // path(ref)

   main:

   contigs_meta
   .filter { meta, tax_info, contigs -> contigs.size() > 0 } 
   .set { contigs_ch }

   REFERENCE_PREP(
    contigs_ch,
    reads_meta,
    ref_ch
   )

   //REFERENCE_PREP.out.ref
   //.join(REFERENCE_PREP.out.contigs, by: [0, 1])
   //.view {"${it}"}

   MUMMER(
    REFERENCE_PREP.out.ref
    .join(REFERENCE_PREP.out.contigs, by: [0, 1])
   )

   FILTER_AND_GLUE_CONTIGS(
    MUMMER.out.delta_tile
   )

   GAPFILL_WITH_READS(
    FILTER_AND_GLUE_CONTIGS
    .out.intermediate_scaffold
    .join(REFERENCE_PREP.out.reads, by: [0, 1])
   )

   GAPFILL_WITH_REF(
    GAPFILL_WITH_READS.out.gapfilled_scaffold
   )

   emit: 
   contigs  = GAPFILL_WITH_REF.out.prep_scaffold
   reads    = REFERENCE_PREP.out.reads 

}
