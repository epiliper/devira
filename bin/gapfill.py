#!/usr/bin/env python3

import subprocess 
import os
import Bio.SeqIO
import argparse
import time
import logging

log = logging.getLogger("gapfill")
log.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
log.addHandler(console_handler)

def run_gap2seq(in_scaffolds, out_scaffolds, reads, kmer_thres, kmer_size, rand_seed):

    cmd = [
            "Gap2Seq",
            "--scaffolds",
            in_scaffolds,
            "--filled",
            out_scaffolds,
            "--reads",
            reads,
            "-k",
            kmer_size,
            "--solid",
            kmer_thres,
            "--randseed",
            rand_seed
            ]

    log.debug(cmd)
    with open("gap2seq.log", "a") as outf:
        subprocess.run(cmd, stdout=outf)

def gapfill(in_scaffold, in_fastq, out_scaffold, solid_kmer_thresholds=[3], kmer_sizes=(90, 80, 70, 60, 50, 40, 31),
                min_gap_to_close=4, time_soft_limit_minutes=60.0, random_seed=0):
        """Try to fill the gaps in the given scaffold, using the reads.

        Inputs:
            in_scaffold: a FASTA file containing the scaffold.  Each FASTA record corresponds to one
                segment (for multi-segment genomes).  Contigs within each segment are
                separated by Ns.  The exact number of Ns between contigs does not matter, as the length of the gap is one 
                of the things determined by the gap-filling tool.  (But see `min_gap_to_close`).
            in_bam: reads to use for filling the gaps.  Only paired-end reads from the bam file are used, any unpaired reads
                are ignored.
           
        Outputs:
            out_scaffold: the input scaffold, with some of the gaps between contigs possibly filled.

        Params:
            solid_kmer_thresholds: kmers must appear at least this many times in the reads to be considered solid.
                We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            kmer_sizes: kmer sizes to use.  We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            min_gap_to_close: stop gap-closing if all gaps are no longer than this many Ns
            gap2seq_opts: extra command-line flags to pass to Gap2Seq
            mem_limit_gb: max memory to use, in gigabytes
            threads: number of threads to use; None means use all available cores.
            time_soft_limit_minutes: stop trying to close more gaps after this many minutes (currently this is a soft/advisory limit)
            random_seed: random seed for choosing random paths (0 to use current time)
        
        """
        stop_time = time.time() + 60* time_soft_limit_minutes
        prev_scaffold = in_scaffold
        print(prev_scaffold)

        for kmer_thres in solid_kmer_thresholds:
            print(kmer_thres)
            for kmer_size in kmer_sizes:

                if not any('N'*min_gap_to_close in str(rec.seq) for rec in Bio.SeqIO.parse(prev_scaffold, 'fasta')):
                    log.info('no gaps left, quittting gap2seq early')
                    break
                if time.time() > stop_time:
                    log.info('Time limit for gap closing reached')
                    break

                filled_scaffold = 'gap2seq-filled.s{}.k{}.fasta'.format(kmer_thres, kmer_size)
                run_gap2seq(in_scaffold, filled_scaffold, in_fastq, str(kmer_thres), str(kmer_size), str(random_seed))

                prev_scaffold = filled_scaffold

        Bio.SeqIO.convert(prev_scaffold, 'fasta', out_scaffold, 'fasta')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_scaffold", help = "The scaffolds you want gapfilled, with Ns in between each contig")
    parser.add_argument("--in_reads", help = "reads to use for gapfilling. Should be a string in the form of '<r1_file>, <r2_file>' if using paired_end data")
    parser.add_argument("--out_scaffold", help = "The output file for the gap-filled scaffold")
    args = parser.parse_args()
    gapfill(args.in_scaffold, args.in_reads, args.out_scaffold)
