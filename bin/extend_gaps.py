#!/usr/bin/env python3

import argparse
# from Bio import SeqIO, pairwise2
from Bio import SeqIO, Align
from Bio.Seq import Seq
import os
import concurrent.futures

def get_seq_by_header(fastas, header):
    for f in fastas:
        if f.id == header:
            return str(f.seq)

    print(f"No fasta found matching header {header}")
    return

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2.0
aligner.mismatch_score =  -1.0
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.5

def align_and_get_overlap(ref_record, query_record):
    ref_id = ref_record.id
    ref_seq = ref_record.seq

    query_id = query_record.id
    query_seq = query_record.seq

    overhang_5prime = 0
    overhang_3prime = 0

    ## this is what this would look like with the old alignment API
    # alignments = pairwise2.align.globalms(
    #     ref_seq, 
    #     query_seq,
    #     2,     # match score
    #     -1,    # mismatch pen
    #     -1.5,    # gap open pen
    #     -0.5,  # gap extend pen
    #     one_alignment_only=True
    # )

    try:
        alignments = aligner.align(ref_seq, query_seq)
        if not alignments:
            return {"query_id": query_id, "overhang_5_prime": overhang_5prime, "overhang_3_prime": overhang_3prime}

    except OverflowError:
        ## crappy alignment
        print("selected reference alignment to contigs is poor, not extending ends.")
        return {"query_id": query_id, "overhang_5_prime": overhang_5prime, "overhang_3_prime": overhang_3prime}

    ## get best alignment
    alignment = alignments[0]

    alignment = alignment.format("fasta")
    ## get just the sequences, not headers or newlines
    alignment = [line for line in alignment.split("\n") if (not line.startswith(">") and line != "")]
    assert len(alignment) == 2

    aligned_ref = alignment[0]
    aligned_query = alignment[1]

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

    return {"query_id": query_id, "overhang_5_prime": overhang_5prime, "overhang_3_prime": overhang_3prime}

def find_max_overhangs(reference_file, query_file, threads):
    """
    find which contigs extend the most past the 5' and 3' ends, respectively, of a reference sequence

    takes in:
    - a file of a single reference fasta (the first fasta will be picked if a multifasta is provided)
    - a multifasta of sequences (recommended to be contigs)
    
    """
    # Load reference sequences
    reference_records = list(SeqIO.parse(reference_file, "fasta"))
    query_records = list(SeqIO.parse(query_file, "fasta"))
    
    if not reference_records:
        raise ValueError("No reference sequences found in the provided file")
    
    # Use the first reference sequence
    ref_record = reference_records[0]
    ref_seq = str(ref_record.seq)
    ref_id = ref_record.id
    
    # a list of [query_id, overhang_5prime, overhang_3prime]
    results = []
    
    # map all contigs to check for possible reference/scaffold extensions
    with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
        futures = []
        for query_record in query_records:
            futures.append(executor.submit(align_and_get_overlap, ref_record, query_record))
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    # sort to contigs with best 5' and 3' overhangs
    results.sort(key=lambda x: x["overhang_5_prime"], reverse=True)
    max_5prime_query = results[0]

    results.sort(key=lambda x: x["overhang_3_prime"], reverse=True)
    max_3prime_query = results[0]
    
    return {
        "reference_id": ref_id,
        "max_5prime_overhang": {
            "query_id": max_5prime_query["query_id"] if max_5prime_query["overhang_5_prime"] > 0 else None,
            "length": max_5prime_query["overhang_5_prime"],
            "seq": get_seq_by_header(query_records, max_5prime_query["query_id"])
        },
        "max_3prime_overhang": {
            "query_id": max_3prime_query["query_id"] if max_3prime_query["overhang_3_prime"] > 0 else None,
            "length": max_3prime_query["overhang_3_prime"],
            "seq": get_seq_by_header(query_records, max_3prime_query["query_id"])
        },
        "all_results": results
    }

def extend_sequence(ref_file, best_5prime, best_3prime, outfile):
        ref_records = list(SeqIO.parse(ref_file, "fasta"))
        ref_id = ref_records[0].id
        ref_seq = new_ref = list(ref_records[0].seq) ## convert to char array

        if best_5prime["query_id"] is not None:
            contig_id = best_5prime["query_id"]
            left_overhang = best_5prime["length"]
            new_ref = list(best_5prime["seq"])[0:left_overhang + 1] + new_ref
            print(f"extended {ref_id} 5' with {left_overhang} bases from contig {contig_id}")

        if best_3prime["query_id"] is not None:
            contig_id = best_3prime["query_id"]
            right_overhang = best_5prime["length"]
            new_ref = list(best_3prime["seq"])[0:right_overhang + 1] + new_ref
            print(f"extended {ref_id} 3' with {right_overhang} bases from contig {contig_id}")


        print(f"length before extension: {len(ref_seq)}; length after extensions: {len(new_ref)}")
        assert len(new_ref) >= len(ref_seq), "Output sequence is shorter than supplied sequence"
        new_ref = "".join(new_ref)

        with open(outfile, "w") as outf:
            outf.write(f">{ref_id}\n")
            outf.write(f"{new_ref}\n")

def main():
    parser = argparse.ArgumentParser(description="Find which query sequence extends past the reference the most")
    parser.add_argument("-r", "--reference", required=True, help="scaffold fasta")
    parser.add_argument("-q", "--query", required=True, help="contig fasta")
    parser.add_argument("-o", "--output", required=True, help="output for extended fasta")
    parser.add_argument("-t", "--threads", required=False, type = int, default=os.cpu_count(), help="number of threads to use")
    
    args = parser.parse_args()
    
    results = find_max_overhangs(args.reference, args.query, args.threads) ## find max overhangs
    extended_ref = extend_sequence(args.reference, results["max_5prime_overhang"], results["max_3prime_overhang"], args.output) ## extend scaffold and write to output

    # debug stuff
    # print(f"Reference: {results['reference_id']}")
    # print("\nMaximum overhangs found:")
    # print(f"5' (left) overhang: {results['max_5prime_overhang']['length']} bp from query {results['max_5prime_overhang']['query_id']}")
    # print(f"3' (right) overhang: {results['max_3prime_overhang']['length']} bp from query {results['max_3prime_overhang']['query_id']}")
    
    # print("\nTop 5 queries by total overhang:")
    # for i, result in enumerate(results['all_results'][:5], 1):
    #     print(f"{i}. {result['query_id']}: {result['total_overhang']} bp total ({result['5prime_overhang']} bp 5', {result['3prime_overhang']} bp 3')")
            
if __name__ == "__main__":
    exit(main())
