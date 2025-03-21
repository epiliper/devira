include { MEGAHIT                   } from '../modules/megahit'
include { METASPADES                } from '../modules/metaspades'
include { FILTER_AND_GLUE_CONTIGS   } from '../modules/filter_and_glue_contigs'
include { EXTEND_SCAFFOLDS          } from '../modules/extend_scaffolds'
include { GAPFILL_WITH_READS        } from '../modules/gapfill_with_reads'
include { GAPFILL_WITH_REF          } from '../modules/gapfill_with_ref'

include { REFERENCE_PREP            } from './reference_prep'

workflow CONTIG_GEN {

   take: 
   reads_meta       // tuple val(meta), val(tax_info), path(reads)
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

   contigs
   .filter { meta, tax_info, contigs -> contigs.size() > 0 } 
   .set { contigs_ch }

   REFERENCE_PREP(
    contigs_ch,
    reads_meta,
    ref_ch
   )

   FILTER_AND_GLUE_CONTIGS(
    REFERENCE_PREP.out.ref
    .join(REFERENCE_PREP.out.contigs, by: [0, 1])
   )

   EXTEND_SCAFFOLDS(
    FILTER_AND_GLUE_CONTIGS
    .out.intermediate_scaffold
    .join(REFERENCE_PREP.out.contigs, by: [0, 1])
   )

   GAPFILL_WITH_READS(
    EXTEND_SCAFFOLDS
    .out.extended_scaffold
    .join(REFERENCE_PREP.out.reads, by: [0, 1])
   )

   emit: 
   contigs      = GAPFILL_WITH_READS.out.gapfilled_scaffold
   reads        = REFERENCE_PREP.out.reads 
   contig_stats = FILTER_AND_GLUE_CONTIGS.out.contig_stats

}
