#!/usr/bin/env python3

import csv
import argparse
from collections import Counter
import subprocess

OUT_COLS = [
    "ref_acc", 
    "ref_tag", 
    "ref_desc", 
    "Query_name", 
    "ANI", 
    "Align_fraction_ref", 
    "Total_bases_covered", 
    ]

def prepare_row(row):
    """
    add reference accession number, tag, and description to the row, and abridge the row to include just the cols in OUT_COLS. 
    Designed to work with fasta reference headers formatted exactly like those in the default ADAR database:
    <acc>SPACE<tag>SPACE<rest of header(AKA desc)>
    """
    s = row["Ref_name"].split(' ')
    acc, tag, desc = s[0], s[1], " ".join(s[2:])

    flu_detected = False

    if "Influenza" in desc:
        flu_detected = True

    new_row = {}
    for col, val in zip(OUT_COLS[:3], [acc, tag, desc]):
        new_row[col] = val

    for col in OUT_COLS[3:]:
        if col not in row:
            new_row[col] = "FORCE-SELECTED"
        else:
            new_row[col] = row[col]

    return flu_detected, new_row

def select_refs(sample_identifier, dist_tsv, ani_thres, align_ref_thres, reference_fasta):
    """
    Parses the TSV output from skani dist to get reference genomes suitable for downstream consensus assembly.
    For every unique reference tag, select the reference with highest % ref covered and ANI, in descending order of priority. Write all qualifiers to TSV.
    If no reference exists within the supplied thresholds, create a failed assembly TSV with the closest one.
    """
    has_flu = False

    with open(dist_tsv, "r") as inf:
        lines = inf.readlines()

    inheader = lines[0].replace("\n", "").split("\t")
    data = list(csv.DictReader(lines[1:], delimiter = "\t", fieldnames=inheader))

    if not data:
        print(f"No references were detected PERIOD (below or above SKANI thresholds) for sample {sample_identifier}")
        return


    # sort SKANI db hits first by fraction of ref genome covered, then by ANI
    data.sort(key = lambda row: (
        float(row["ANI"]),
        float(row["Total_bases_covered"]),
        float(row["Ref_50_ctg_len"])
        ), reverse = True)

    out_rows = {}
    for row in data:
        s = row["Ref_name"].split(' ')
        acc, tag, desc = s[0], s[1], s[2:]
        if tag not in out_rows:
            if (float(row["Align_fraction_ref"]) >= align_ref_thres and 
                float(row["ANI"]) >= ani_thres):

                flu_detected, out_row = prepare_row(row)

                if flu_detected:
                    has_flu = True

                out_rows[tag] = out_row

    if not out_rows: ## no refs found passing thresholds
        print(f"No suitable references found for sample {sample_identifier}!")
        out_rows["_"] = prepare_row(data[0])
        report_file = f"{sample_identifier}_failed_assembly.tsv"

    else:
        print(f"Found at least one reference to use for {sample_identifier}.")
        report_file = f"{sample_identifier}_covstats.tsv"

        if has_flu:
            print("FLU DETECTED")
            ## if we have a match to fewer than 8 flu segments, force-select all remaining segments for that subtype
            ## this assumes ref tag is formatted as H*N*_SEGMENT_NAME (for Flu A) or Flu_TYPE_SEGMENT_NAME (for Flus B, C, D)
            ## if fewer than 8 segments for most common subtype, add missing segments, keep the rest

            ## get the most prevalent subtype
            counter = Counter()
            for r in out_rows:
                counter[r[:4]] += 1

            subtype, count = counter.most_common(1)[0]
            if count < 8: # missing segments
                print(f"Attempting for force-select all 8 flu segments for {sample_identifier}")

                cmd = ["grep", subtype, reference_fasta]
                hits = subprocess.check_output(cmd, text = True)
                if hits:
                    hits = hits.split("\n")[:-1] # get rid of blank line
                    print(hits)
                    for header in hits:
                        header = header[1:] # get rid of >
                        tag = header.split(" ")[1]
                        if tag not in out_rows:
                            _, out_row = prepare_row({"Ref_name": header})
                            out_rows[tag] = out_row

    with open(report_file, "w") as outf:
        writer = csv.DictWriter(outf, fieldnames = OUT_COLS, delimiter = "\t")
        writer.writeheader()
        for row in out_rows.values():
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--sample", help = "A descriptor for the input file, used for naming output", required = True)
    _ = parser.add_argument("--dist_report", help = "The output TSV from skani dist", required = True)
    _ = parser.add_argument("--ani_thres", type = float, help = "Threshold for ANI", required = True)
    _ = parser.add_argument("--align_ref_thres", type = float, help = "Threshold for % ref covered by query", required = True)
    _ = parser.add_argument("--reference_fasta", type = str, help = "Reference fasta used for force-selecting needed reference headers")
    args = parser.parse_args()

    select_refs(args.sample, args.dist_report, args.ani_thres, args.align_ref_thres, args.reference_fasta)

