#!/usr/bin/env python3

import pysam
import sys
from Bio import SeqIO
from Bio.Seq import reverse_complement
import argparse
import logging
import gapfill
import random
from enum import Enum

class Extend(Enum):
    LEFT = "5'"
    INTERIOR = "interior"
    RIGHT = "3'"

log = logging.getLogger("scaffold")
log.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
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
    if not query_seqs:
        longest = num_filtered = shortest = n50 = scaf_len = -1

    else:
        lengths = [len(seq) for seq in query_seqs]
        lengths.sort(reverse=True)
        longest = lengths[0]
        shortest = lengths[-1]

        def get_N50(lengths: list[int]):
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
    ref_len: int
    aln_record: pysam.AlignedSegment
    id: str
    length: int
    aln_start: int
    query_start: int
    aln_end: int
    query_end: int

    def has_left_overhang(self) -> int:
        """
        check if the alignment has clipped sequence that precedes the 0th ref coordinate
        if it does, return the number of bases from 0th coordinate to the alignment start
        """
        rec = self.aln_record
        cig = rec.cigartuples
        if not cig: return 0

        if cig[0][0] in [4, 5]:
            clip_len = cig[0][1]
            if self.aln_start - clip_len < 0:
                return self.aln_start

        return 0

    def has_right_overhang(self) -> int:
        """
        check if the alignment has clipped sequence that goes beyond the last ref coordinate
        if it does, return the number of bases from the end of the alignment to the end of the reference
        """
        rec = self.aln_record
        cig = rec.cigartuples
        if not cig or len(cig) == 1: return 0

        if cig[-1][0] in [4, 5]:
            clip_len = cig[-1][1]
            if self.aln_end + clip_len >= self.ref_len:
                return self.ref_len - self.aln_end

        return 0

    def trim_seq_to_aln(self):
        self.seq = self.seq[self.query_start: self.query_end]
        assert len(self.seq) == self.length
                
    def trim_seq_to_ref(self):
        l_delta = self.aln_start
        r_delta = self.aln_end - self.ref_len

        if l_delta < 0:
            l_delta = abs(l_delta)
            self.aln_start += l_delta
            self.length -= l_delta
            self.query_start += l_delta
            assert self.aln_start == 0

        if r_delta > 0:
            self.aln_end -= r_delta
            self.length -= r_delta
            self.query_end -= r_delta
            assert(self.aln_end == self.ref_len)

        self.seq = self.seq[self.query_start: self.query_end]

    def extend_clip(self, ext_choice: Extend, limit: int):
        rec = self.aln_record
        cig = rec.cigartuples
        did_extension = False
        pad = 0
        if not cig: return

        match ext_choice:
            case Extend.LEFT:
                if cig[0][0] in [4, 5]:
                    pad = min(cig[0][1], limit)
                    self.length += pad
                    self.aln_start -= pad
                    self.query_start -= pad
                    did_extension = True

            case Extend.INTERIOR:
                if len(cig) <= 2:
                    return

                for c in cig[1:-1]:
                    # interior clipping detected
                    if c[0] in [4, 5]:
                        pad = cig[0][1]
                        self.length += pad
                        self.aln_end += pad
                        self.query_end += pad
                        did_extension = True

            case Extend.RIGHT:
                if len(cig) == 1:
                    return

                if len(cig) > 1 and cig[-1][0] in [4, 5]:
                    pad = min(cig[-1][1], limit) # limit end ext
                    self.length += pad
                    self.aln_end += pad
                    self.query_end += pad
                    did_extension = True

        if self.length > len(self.seq):
            log.fatal(f"{self.id}: Error: cigar length is longer than query sequence: {self.length}: {len(self.seq)}")
            exit(1)

        if did_extension:
            log.info(f"{self.id}: extended to length {self.length} from cigar {cig} along coordinates {self.aln_start} - {self.aln_end}; Extension type: {ext_choice.value}")

    def load_seq(self, seq: str, rev: bool):
        if rev:
            self.seq = str(reverse_complement(seq))
        else:
            self.seq = seq

    def __init__(self, rec: pysam.AlignedSegment, seq: str, ref_len: int):
        self.aln_record = rec
        self.ref_len = ref_len
        self.id = rec.query_name or "_NAME_NOT_FOUND_"
        self.length = self.aln_record.query_alignment_length
        self.query_start = self.aln_record.query_alignment_start
        self.query_end = self.aln_record.query_alignment_end

        self.aln_start = self.aln_record.reference_start
        ## calculate reference end if it's not available
        end = self.aln_record.reference_end 
        if not end:
            end = self.aln_record.query_alignment_length + self.aln_record.reference_start
        self.aln_end = end

        self.load_seq(seq, rec.is_reverse)

        # if the sequence overhangs the ref on the left, then limit extension to either the entire overhang, or 1/10th ref len, whichever is longest
        # if no overhang, then don't limit extension
        l_overhang_limit = max(self.has_left_overhang(), int(self.ref_len * 0.10))

        # if the sequence overhangs the ref on the right, then limit extension to either the entire overhang, or 1/10th ref len, whichever is longest
        # if no overhang, then don't limit extension
        r_overhang_limit = max(self.has_right_overhang(), int(self.ref_len * 0.10))

        self.extend_clip(Extend.LEFT, l_overhang_limit)
        self.extend_clip(Extend.INTERIOR, sys.maxsize)
        self.extend_clip(Extend.RIGHT, r_overhang_limit)

        self.trim_seq_to_aln()

def load_alignments(file: str, fastas: dict[str, str], min_seq_length: int) -> tuple[int, list[AlignedSequence]]:
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
            log.info(f"Skipping loading alignment {rec.query_name}; mapped: {rec.is_mapped}; supplementary aln: {rec.is_supplementary}")
            continue

        try:

            seq = fastas[rec.query_name] 
            if len(seq) < min_seq_length:
                log.info(f"Skipping loading alignment for {rec.query_name}; contig sequence is too short: {len(seq)} / {min_seq_length} ")
                continue

        except KeyError:
            log.info(f"Sequence for alignment not found: {rec.query_name}")
            exit(1)

        records.append(AlignedSequence(rec, seq, ref_len))

    return (ref_len, records)

def load_fastas(file: str) -> tuple[int, dict[str, str]]:
    """
    creates a dict of <fasta_id>:<fasta_sequence> mappings for each entry in a given fasta file

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
    if not recs:
        return None

    recs.sort(key = lambda x: (x.length, x.aln_record.mapping_quality), reverse = True)
    ret: AlignedSequence | None = None
    r = recs[0]

    ## this contig must span from the start of (or before) the reference, until the end or after
    if r.length >= ref_len and r.aln_start <= 0:
        ret = r

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
    if not recs:
        return "".join(["N"] * ref_len)

    supercontig = check_if_longer_than_ref(recs, ref_len)
    if supercontig:
        log.info(f"{supercontig.id}:Found supercontig for reference, from coords {supercontig.aln_start} - {supercontig.aln_end}")
        seq = supercontig.seq[0:supercontig.length]
        return seq

    seq = []
    recs.sort(key = lambda x: x.aln_start)
    start, end = recs[0].aln_start, 0
    num_ns = 0
    
    for r in recs:
        if r.aln_start > end:
            # a gap exists between this contig and either 1) the 0th reference coord or 2) the previous contig.
            add = ["N"] * (r.aln_start - end) + list(r.seq[0: r.length])
            num_ns += (r.aln_start - end)

        elif r.aln_start < end and end != 0:
            # part of this contig has already been covered by the previous one.
            add = list(r.seq)[end - r.aln_start:]

        else:
            # this contig starts immediately after the previous contig (or is the first one)
            add = list(r.seq)

        seq += add
        delta = len(add)

        log.info(f"{r.id}: glued {delta} bases; Sequence has {num_ns} Ns...")

        # make sure end only keeps track of bases added after ref coord 0; it shouldn't count bases added before
        delta += min(0, r.aln_start)

        end += delta

    if pad_ends:
        if start > 0:
            seq = ["N"] * start + seq

        if end < ref_len:
            seq = seq + ["N"] * (ref_len - end)

    return "".join(seq)


def pad_with_ns(seq: str, l_pad: int, r_pad: int):
    log.info(f"extending sequence; 5' with {l_pad} Ns; 3' with {r_pad} Ns")
    return "".join(["N"] * l_pad + list(seq) + ["N"] * r_pad)


def main(align_file: str, query_fasta: str, outfile: str, prefix: str, reads: str, report_file: str, min_seq_length: int):
    log.info(f"Filtering contigs used for tiling to length >= {min_seq_length} bp")
    num_fastas, fasta_dict = load_fastas(query_fasta)
    ref_len, aligns = load_alignments(align_file, fasta_dict, min_seq_length)

    for a in aligns:
        log.debug(a.seq, a.length, a.aln_start, a.aln_end, a.aln_record.query_name)

    seq = glue_alns_across_ref(aligns, ref_len, False)

    ## if we have a significant number of continuous Ns, then try to gapfill with Gap2Seq
    seq = seq.lstrip("N").rstrip("N") # strip leading and trailing Ns, which otherwise make gap2seq crash
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
    _ = parser.add_argument("-l", "--min_contig_length", type = int, required = False, default = 300, help = "minimum length of raw contig sequence for contig to be used in scaffolding")
    args, _ = parser.parse_known_args()

    main(args.alignment, args.query, args.output, args.prefix, args.reads, args.contig_report, args.min_contig_length)
