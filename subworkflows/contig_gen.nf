include { MEGAHIT                   } from '../modules/megahit'
include { METASPADES                } from '../modules/metaspades'
include { MUMMER                    } from '../modules/mummer'
include { FILTER_AND_GLUE_CONTIGS   } from '../modules/filter_and_glue_contigs'
include { GAPFILL_WITH_READS        } from '../modules/gapfill_with_reads'
include { GAPFILL_WITH_REF          } from '../modules/gapfill_with_ref'

include { REFERENCE_PREP            } from './reference_prep'

workflow CONTIG_GEN {

   take: 
   reads_meta       // tuple val(meta), path(reads)
   contig_method    // val(contig_method)
   ref_ch           // path(ref)

   main:

   if (contig_method == "megahit") {
       MEGAHIT(reads_meta)
    contigs = MEGAHIT.out.contigs

   } 
   else if (contig_method == "metaspades") {
       METASPADES(reads_meta)
    contigs = METASPADES.out.contigs
   }

   else {
    // this should be checked beforehand in main.nf, but leaving it here for completion's sake
    throw new IllegalArgumentException("Invalid contig method: ${contig_method}. Choose from 'metaspades' or 'megahit'")
   }

   REFERENCE_PREP(
    contigs,
    reads_meta,
    ref_ch
   )

   MUMMER(
    REFERENCE_PREP.out.ref
    .join(REFERENCE_PREP.out.contigs)
   )

   FILTER_AND_GLUE_CONTIGS(
    MUMMER.out.delta_tile
   )

   GAPFILL_WITH_READS(
    FILTER_AND_GLUE_CONTIGS
    .out.intermediate_scaffold
    .join(REFERENCE_PREP.out.reads)
   )

   GAPFILL_WITH_REF(
    GAPFILL_WITH_READS.out.gapfilled_scaffold
   )

   emit: 
   contigs  = GAPFILL_WITH_REF.out.prep_scaffold
   reads    = REFERENCE_PREP.out.reads 

}
