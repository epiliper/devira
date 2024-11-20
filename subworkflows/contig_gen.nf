include { MEGAHIT } from '../modules/megahit'
include { METASPADES } from '../modules/metaspades'

workflow CONTIG_GEN {

   take: 
   reads_meta // tuple val(meta), path(reads)
   contig_method // val(contig_method)

   main:

   if (contig_method == "megahit") {

       MEGAHIT(
               reads_meta,
              )

    contigs = MEGAHIT.out.contigs_fasta
    sid = MEGAHIT.out.sample_id

   } 

   else if (contig_method == "metaspades") {

       METASPADES(
               reads_meta,
               )
    
    contigs = METASPADES.out.contigs_fasta
    sid = METASPADES.out.sample_id

   }

   else {

    throw new IllegalArgumentException("Invalid contig method: ${contig_method}. Choose from 'metaspades' or 'megahit'")
   }

   emit: 
   contigs = contigs
   sid     = sid

}
