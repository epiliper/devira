#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def get_filtered_contigs(tile_file):
    contigs = []
    with open(tile_file) as inf:
        for line in inf:
            t = line.split("\t")
            contigs.append(t[12].replace("\n", ""))

    return contigs

def load_contigs(tile_file, contig_file):
    contig_count = 0
    filtered_contig_ids = get_filtered_contigs(tile_file)
    queries = list(SeqIO.parse(contig_file, "fasta"))
    filtered_contigs = [q for q in queries if q.id in filtered_contig_ids]
    return queries, filtered_contigs

def report_contig_stats(query_records):
    if not query_records:
        return -1, -1, -1, -1

    lengths = [len(record.seq) for record in query_records]
    lengths.sort(reverse=True)
    longest = lengths[0]
    shortest = lengths[-1]

    def get_N50(query_lengths):
        total_length = sum(lengths)

        cumu_length = 0
        for length in lengths:
            cumu_length += length
            if cumu_length >= total_length / 2:
                return length

        return -1

    n50 = get_N50(lengths)
    num_contigs = len(query_records)
    return num_contigs, shortest, longest, n50



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--contigs", help="(multi)fasta of contigs")
    parser.add_argument("-t", "--tile", help="Mummer tile file")
    parser.add_argument("-o", "--out", help="name of output tsv")
    args = parser.parse_args()

    total, filtered = load_contigs(args.tile, args.contigs)
    num_total = len(total)
    num_filtered, shortest, longest, n50 = report_contig_stats(filtered)

    if not args.out.endswith("tsv"): args.out += ".tsv"

    with open(args.out, "w") as outf:
        outf.write("num_total_contigs\tnum_filtered_contigs\tN50\tshortest_contig_len\tlongest_contig_len\n")
        outf.write(f"{num_total}\t{num_filtered}\t{n50}\t{shortest}\t{longest}")
