#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Align
import os
import concurrent.futures
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

class AlnResult:
    id: str | None
    overhang5: int
    overhang3: int
    seq: str

    def __init__(self, id:str | None, overhang_5:int, overhang_3:int, seq:str, reverse = False):
        self.id = id
        self.overhang5 = overhang_5
        self.overhang3 = overhang_3
        self.seq = str(reverse_complement(seq)) if reverse else seq

    def get_5prime_overhang(self, n_fill: bool = False):
        if n_fill:
            over = ["N" for _ in range(0, self.overhang5)]
        else:
            over = list(self.seq)[0:self.overhang5]

        return len(over), over

    def get_3prime_overhang(self, n_fill: bool = False):
        overhang_3_start_pos = len(self.seq) - 1 - self.overhang3
        if n_fill:
            over = ["N" for _ in range(overhang_3_start_pos, len(self.seq))]
        else:
            over = list(self.seq)[overhang_3_start_pos: ]
        return len(over), over

    def get(self, attr: str) -> int:
        return getattr(self, attr)

class AlnSummary:
    ref_id: str | None
    max_5_overhang: AlnResult | None
    max_3_overhang: AlnResult | None
    all_alns: list[AlnResult]
    
    def __init__(self, ref_id: str | None, all_alns: list[AlnResult], thres: float):
        self.ref_id = ref_id
        self.all_alns = all_alns
        self.max_5_overhang = self.get_max_below_thres("overhang5", thres)
        self.max_3_overhang = self.get_max_below_thres("overhang3", thres)

    def get_max_below_thres(self, attr: str, thres: float) -> AlnResult | None:
        self.all_alns.sort(key = lambda x: x.get(attr), reverse=True)
        found = False
        q = None
        for q in self.all_alns:
            if getattr(q, attr) <= thres:
                found = True
                break

        if not found or not q:
            return None

        return q

aligner = Align.PairwiseAligner(scoring = "blastn")
aligner.mode = 'global'

def load_fasta(fasta_file: str, load_first: bool) -> list[SeqRecord]:
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) > 1 and load_first:
        print(f"Warning: loading first record of multiple records in file {fasta_file}")

    if load_first:
        return [records[0]]

    return records


def align_and_get_overlap(ref_record: SeqRecord, query_record: SeqRecord):
    dud = AlnResult("", -1, -1, "") ## return in case no alignments found
    ref_seq = ref_record.seq

    query_id = query_record.id
    query_seq = query_record.seq

    overhang_5prime = 0
    overhang_3prime = 0
    reverse = False

    try:
        alignments = aligner.align(ref_seq, query_seq, strand = "+")
        if not alignments:
            return dud

    except OverflowError:
        ## crappy alignment
        print(f"{query_id}: reference alignment is poor with forward strand. Trying reverse...")

    try:
        alignments = aligner.align(ref_seq, query_seq, strand = "-")
        ## if we keep this alignment, we need to remember to reverse-complement the query sequence 
        ## for downstream extension
        reverse = True 
        if not alignments:
            return dud

    except OverflowError:
        print(f"{query_id}: no suitable alignments found for either strand. Skipping...")
        return dud

    ## get best alignment
    alignment = alignments[0]

    if alignment.score <= 500:
        print(f"{query_id}: no suitable alignments found for either strand. Skipping...")
        return dud

    print(f"Proceeding with alignment of score {alignment.score}")

    alignment = alignment.format("fasta")
    ## get just the sequences, not headers or newlines
    alignment = [line for line in alignment.split("\n") if (not line.startswith(">") and line != "")]
    assert len(alignment) == 2

    aligned_ref: str = alignment[0]
    aligned_query: str = alignment[1]

    # get 5' overhang
    overhang_5prime = 0
    for i in range(len(aligned_ref)):
        if aligned_ref[i] == '-':
            if aligned_query[i] != '-':
                overhang_5prime += 1
        else:
            break
    
    # get 3' overhang
    overhang_3prime = 0
    for i in range(len(aligned_ref) - 1, -1, -1):
        if aligned_ref[i] == '-':
            if aligned_query[i] != '-':
                overhang_3prime += 1
        else:
            break

    return AlnResult(query_id, overhang_5prime, overhang_3prime, query_seq, reverse)

def find_max_overhangs(rec_to_extend: SeqRecord, query_recs: list[SeqRecord], threads: int) -> AlnSummary:
    """
    find which contigs extend the most past the 5' and 3' ends, respectively, of a reference sequence

    takes in:
    - a file of a single reference fasta (the first fasta will be picked if a multifasta is provided)
    - a multifasta of sequences (recommended to be contigs)
    
    """
    # Use the first reference sequence
    original_seq = str(rec_to_extend.seq)
    input_id = rec_to_extend.id
    
    # a list of [query_id, overhang_5prime, overhang_3prime]
    results = []
    
    # map all contigs to check for possible reference/scaffold extensions
    with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
        futures = []
        for query_record in query_recs:
            futures.append(executor.submit(align_and_get_overlap, rec_to_extend, query_record))
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())


    aln_summary = AlnSummary(input_id, results, len(original_seq) / 3)
    return aln_summary
   
def extend_seq_w_overhangs(rec_to_extend: SeqRecord, results: AlnSummary, n_fill: bool = False) -> SeqRecord:
        original_seq = new_seq = list(rec_to_extend.seq) ## convert to char array
        ref_id = rec_to_extend.id

        if results.max_5_overhang:
            best_5prime = results.max_5_overhang
            ext_count, ext = best_5prime.get_5prime_overhang(n_fill)
            new_seq = ext + new_seq
            print(f"extended {ref_id} 5' with {ext_count} bases from contig {best_5prime.id}")

        if results.max_3_overhang:
            best_3prime = results.max_3_overhang
            ext_count, ext = best_3prime.get_3prime_overhang(n_fill)
            new_seq = new_seq + ext
            print(f"extended {ref_id} 3' with {ext_count} bases from contig {best_3prime.id}")


        print(f"length before extension: {len(original_seq)}; length after extensions: {len(new_seq)}")
        assert len(new_seq) >= len(original_seq), "Output sequence is shorter than supplied sequence"
        new_seq = "".join(new_seq)

        new_rec = SeqRecord(id=ref_id, seq=Seq(data = new_seq, length=len(new_seq)))

        return new_rec


def main():
    parser = argparse.ArgumentParser(description="Find which query sequence extends past the reference the most")
    _ = parser.add_argument("-s", "--scaffold", required=True, type=str, help="scaffold fasta")
    _ = parser.add_argument("-r", "--reference", required=True, type=str, help="reference fasta")
    _ = parser.add_argument("-q", "--query", required=True, type=str, help="contig fasta")
    _ = parser.add_argument("-o", "--output", required=True, type=str, help="output for extended fasta")
    _ = parser.add_argument("-t", "--threads", required=False, type=int, default=os.cpu_count(), help="number of threads to use")
    args = parser.parse_args()

    ref = load_fasta(args.reference, load_first=True)[0]
    scaffold = load_fasta(args.scaffold, load_first=True)[0]
    contigs = load_fasta(args.query, load_first=False)

    ## attempt to extend scaffold with contigs
    ## if scaffold is still shorter than reference, extend scaffold ends with Ns up to reference ends
    results: AlnSummary = find_max_overhangs(scaffold, contigs, args.threads)
    extended_record = extend_seq_w_overhangs(scaffold, results)

    if len(extended_record.seq) < len(ref.seq):
        print("assembly is still shorter than reference. Attempting to extend with reference. Using Ns...")
        results: AlnSummary = find_max_overhangs(extended_record, [ref], args.threads)
        extended_record = extend_seq_w_overhangs(extended_record, results)
    else:
        print("assembly longer than reference. Won't extend with reference.")

    with open(args.output, "w") as outf:
        outf.write(f">{extended_record.id}\n")
        outf.write(f"{extended_record.seq}\n")


if __name__ == "__main__":
    exit(main())
