#!/usr/bin/env python3

import pysam
from sys import exit, version
from Bio import SeqIO
from Bio.Seq import reverse_complement
import argparse
import logging
import gapfill
import random

log = logging.getLogger("scaffold")
log.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
log.addHandler(console_handler)

def report_contig_stats(query_seqs: list[str], num_fastas: int, scaf_len: int, output_file: str):
    """
    Reports N50, number of total and filtered fastas, and scaffold len to an output TSV file.

    Args:
        query_seqs: a list of seqs in string type
        num_fastas: the number of seqs prior to any filtering
        scaf_len: length of scaffold assembled with the seqs
        output file: output tsv to write to.

    Returns:
        None
    """
    lengths = [len(seq) for seq in query_seqs]
    lengths.sort(reverse=True)
    longest = lengths[0]
    shortest = lengths[-1]

    def get_N50(query_lengths: list[int]):
        """
        Get N50 from a list of contig lengths sorted in descending order
        """
        total_length = sum(lengths)

        cumu_length = 0
        for length in lengths:
            cumu_length += length
            if cumu_length >= total_length / 2:
                return length

        return -1

    n50 = get_N50(lengths)
    num_filtered = len(query_seqs)

    log.info(f"Writing contig stats to report file: {output_file}")
    with open(output_file, "w") as outf:
        _ = outf.write("num_total_contigs\tnum_filtered_contigs\tN50\tshortest_contig_len\tlongest_contig_len\tscaffold_length\n" + 
                       f"{num_fastas}\t{num_filtered}\t{n50}\t{shortest}\t{longest}\t{scaf_len}\r")

class AlignedSequence:
    seq: str
    aln_record: pysam.AlignedSegment
    length: int
    start: int
    end: int

    def extend_with_clipping(self):
        rec = self.aln_record
        length = rec.query_alignment_length
        cig = rec.cigartuples
        self.start = self.aln_record.reference_start - 1 or 0
        self.end = self.aln_record.reference_end - 1 or 0

        if not cig: return

        ## if the alignment contains soft- (4) or hard-clipping (5), adjust the length to include clipped bases.
        ## We're doing this because we want as little reference bias as possible, and don't want to trim contig sequence yet.
        if cig[0][0] in [4, 5]:
            length += cig[0][1]
            self.start -= cig[0][1]

        for c in cig[1:]:
            if c[0] in [4, 5]:
                length += c[1]
                log.debug(c[0], c[1])
                self.end += c[1]

        self.length = length
        if self.length > len(self.seq):
            log.fatal(f"Error: cigar length is longer than query sequence: {self.length}: {len(self.seq)}")
            exit(1)

    def load_seq(self, seq: str, rev: bool):
        if rev:
            self.seq = str(reverse_complement(seq))
        else:
            self.seq = seq

    def __init__(self, rec: pysam.AlignedSegment, seq: str):
        self.aln_record = rec
        self.load_seq(seq, rec.is_reverse)
        self.extend_with_clipping()


def load_alignments(file: str, fastas: dict[str, str]) -> tuple[int, list[AlignedSequence]]:
    """
    Parses a bam file, loads all mapped and primary alignments into [AlignedSequence] with corresponding query sequence

    Args:
        file: alignment file (sam or bam)
        fastas: a dictionary mapping fasta headers to fasta sequences; used to associate sequence with each sam record

    Returns:
        int: number of mapped, primary alignments
        list[AlignedSequence]: the loaded alignment records

    Raises:
        KeyError if bam record id is not present in fasta dict
    """
    records: list[AlignedSequence]= []
    if file.endswith(".bam"):
        read_mode = "rb"
    elif file.endswith(".sam"):
        read_mode = "r"
    else:
        exit("Invalid alignment file extension: should be .bam or .sam")

    aln_file = pysam.AlignmentFile(file, read_mode)

    if len(aln_file.references) != 1:
        log.fatal("Alignment file should only have one reference. Exiting...")
        exit(1)

    ref_len = aln_file.get_reference_length(aln_file.references[0])

    for rec in aln_file.fetch():
        if not rec.is_mapped or rec.is_supplementary:
            continue
        try:
            seq = fastas[rec.query_name] 
        except KeyError:
            log.fatal(f"sequence for {rec.query_name} not found in bam file. Exiting...")
            exit(1)

        records.append(AlignedSequence(rec, seq))

    return (ref_len, records)

def load_fastas(file: str) -> tuple[int, dict[str, str]]:
    """
    creates a <fasta_id>:<fasta_sequence> dict for each entry in a given fasta file

    Args: 
        file: fasta file

    Returns:
        int: number of fasta records loaded
        dict[str, str]: dictionary mapping fasta header to fasta seq
    """
    fasta_dict = { record.id : str(record.seq) for record in SeqIO.parse(file, "fasta") }
    num_fastas = len(fasta_dict)
    return num_fastas, fasta_dict

def check_if_longer_than_ref(recs: list[AlignedSequence], ref_len: int) -> AlignedSequence | None:
    """
    a simple function to check if any alignment records cover or extend over the entirety of a reference sequence

    Args:
        recs: a list of candidate queries to check
        ref_len: the length of the reference in question

    Returns:
        AlignedSequence | None: a contig that covers all or more of the reference, if it exists, or else None.

    Note that the "length" here is intended to be calculated from # of reference bases covered; this function shouldn't pick partial alignments with massive overhangs.
    """
    recs.sort(key = lambda x: (x.length, x.aln_record.mapping_quality), reverse = True)
    ret: AlignedSequence | None = None
    r = recs[0]

    if r.length > ref_len:
        return ret

def glue_alns_across_ref(recs: list[AlignedSequence], ref_len: int, pad_ends: bool = False) -> str:
    """
    A simple function to orient contigs along a reference, with preference to contig insertions/deletions compared to the reference.

    Args:
        recs: contigs to scaffold
        ref_len: length of reference
        pad_ends: if true, pad end gaps of the scaffold with Ns

    Returns:
        str: the sequence of the scaffold
    """
    supercontig = check_if_longer_than_ref(recs, ref_len)
    if supercontig:
        return supercontig.seq

    seq = []
    recs.sort(key = lambda x: x.start)
    start, end = recs[0].start, -1 
    prev_len = 0
    num_ns = 0
    
    for r in recs:
        if r.start > end and end != -1:
            # a gap exists between this contig and the previous one.
            seq += ["N"] * (r.start - end)
            seq += list(r.seq)[0: r.length]
            num_ns += (r.start - end)

        elif r.start < end and end != -1:
            # part of this contig has already been covered by the previous one.
            seq += list(r.seq)[end - r.start:]

        else:
            seq += list(r.seq)

        delta = len(seq) - prev_len
        prev_len = end = prev_len + delta

        log.info(f"glued {delta} bases; Sequence has {num_ns} Ns...")

    if pad_ends:
        if start > 0:
            seq = ["N"] * start + seq

        if end < ref_len:
            seq = seq + ["N"] * (ref_len - end)

    return "".join(seq)


def pad_with_ns(seq: str, l_pad: int, r_pad: int):
    log.info(f"extending sequence; 5' with {l_pad} Ns; 3' with {r_pad} Ns")
    return "".join(["N"] * l_pad + list(seq) + ["N"] * r_pad)


def main(align_file: str, query_fasta: str, outfile: str, prefix: str, reads: str, report_file: str):
    num_fastas, fasta_dict = load_fastas(query_fasta)
    ref_len, aligns = load_alignments(align_file, fasta_dict)

    for a in aligns:
        log.debug(a.seq, a.length, a.start, a.end, a.aln_record.query_name)

    ## do the scaffolding with just the contigs, don't fill end gaps
    seq = glue_alns_across_ref(aligns, ref_len, False)

    ## if we have a significant number of continuous Ns, then try to gapfill with Gap2Seq
    if "".join(["N"] * 10) in seq:

        temp = f"{random.randint(1_000_000, 9_999_999)}.fa"

        with open(temp, "w") as outf:
            _ = outf.write(f">{prefix}_scaffold\n")
            _ = outf.write(f"{seq}\n")

        gapfill.gapfill(temp, reads, outfile)

        seq = next(SeqIO.parse(outfile, "fasta")).seq

    ## if the gapfilled scaffold is still shorter, then pad it with Ns for downstream mapping/consensus calling
    if len(seq) < ref_len:
        log.info(f"Scaffold is still shorter than reference: {len(seq)}: {ref_len}")
        seq = pad_with_ns(seq, ref_len - len(seq), ref_len - len(seq))

        ## sanity check
        if len(seq) < ref_len: 
            log.fatal(f"Error: scaffold length shorter than reference {len(seq)}, {ref_len}. Aborting...")
            exit(1)

    ## write final sequence to output file.
    with open(outfile, "w") as outf:
        _ = outf.write(f">{prefix}_scaffold\n")
        _ = outf.write(f"{seq}\n")

    report_contig_stats([str(a.seq) for a in aligns], num_fastas, len(seq), report_file)
    log.info(f"Wrote scaffold of length {len(seq)}, made with reference of length {ref_len}, to file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("-a", "--alignment", type = str, help = "alignment sam or bam")
    _ = parser.add_argument("-q", "--query", type = str, help = "query fasta")
    _ = parser.add_argument("-r", "--reads", type = str, help = "reads to use for gapfilling. if mulitple files, should be comma-delimited")
    _ = parser.add_argument("-o", "--output", type = str, help = "name of output fasta file")
    _ = parser.add_argument("-p", "--prefix", type = str, help = "name to use for scaffold")
    _ = parser.add_argument("-c", "--contig_report", type = str, help = "name of contig stats file")
    args, _ = parser.parse_known_args()

    main(args.alignment, args.query, args.output, args.prefix, args.reads, args.contig_report)
