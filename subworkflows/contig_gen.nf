include { MEGAHIT                                   } from '../modules/megahit'
include { METASPADES                                } from '../modules/metaspades'
include { ALIGN_CONTIGS_TO_REF                      } from '../modules/align_contigs_to_ref'
include { FILTER_AND_GLUE_CONTIGS                   } from '../modules/filter_and_glue_contigs'
include { REFERENCE_PREP                            } from './reference_prep'

workflow CONTIG_GEN {

   take: 
   reads_meta       // tuple val(meta), val(tax_info), path(reads)
   contig_method    // val(contig_method)
   ref_ch           // path(ref)
   min_contig_length // val(length)

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

   ALIGN_CONTIGS_TO_REF(
    REFERENCE_PREP.out.contigs.join(
    REFERENCE_PREP.out.ref, by: [0, 1])
   )

   FILTER_AND_GLUE_CONTIGS(
    ALIGN_CONTIGS_TO_REF.out.alignment
    .join(REFERENCE_PREP.out.contigs, by: [0, 1])
    .join(REFERENCE_PREP.out.reads, by: [0, 1]),
    min_contig_length
   )

   emit: 
   contigs      = FILTER_AND_GLUE_CONTIGS.out.intermediate_scaffold
   reads        = REFERENCE_PREP.out.reads 
   contig_stats = FILTER_AND_GLUE_CONTIGS.out.contig_stats

}
