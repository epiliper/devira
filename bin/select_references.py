#!/usr/bin/env python3

import csv
import argparse

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

    new_row = {}
    for col, val in zip(OUT_COLS[:3], [acc, tag, desc]):
        new_row[col] = val

    for col in OUT_COLS[3:]:
        new_row[col] = row[col]

    return new_row

def select_refs(sample_identifier, dist_tsv, ani_thres, align_ref_thres):
    """
    Parses the TSV output from skani dist to get reference genomes suitable for downstream consensus assembly.
    For every unique reference tag, select the reference with highest % ref covered and ANI, in descending order of priority. Write all qualifiers to TSV.
    If no reference exists within the supplied thresholds, create a failed assembly TSV with the closest one.
    """
    with open(dist_tsv, "r") as inf:
        lines = inf.readlines()

    inheader = lines[0].replace("\n", "").split("\t")
    data = list(csv.DictReader(lines[1:], delimiter = "\t", fieldnames=inheader))

    # sort SKANI db hits first by fraction of ref genome covered, then by ANI
    data.sort(key = lambda row: (
        float(row["Total_bases_covered"]),
        float(row["ANI"]),
        float(row["Align_fraction_ref"])
        ), reverse = True)

    out_rows = {}
    for row in data:
        s = row["Ref_name"].split(' ')
        acc, tag, desc = s[0], s[1], s[2:]
        if tag not in out_rows:
            if (float(row["Align_fraction_ref"]) >= align_ref_thres and 
                float(row["ANI"]) >= ani_thres):

                out_row = prepare_row(row)
                out_rows[tag] = out_row

    if not out_rows: ## no refs found passing thresholds
        print(f"No suitable references found for sample {sample_identifier}!")
        out_rows["_"] = prepare_row(data[0])
        report_file = f"{sample_identifier}_failed_assembly.tsv"

    else:
        print(f"Found at least one reference to use for {sample_identifier}.")
        report_file = f"{sample_identifier}_covstats.tsv"

    with open(report_file, "w") as outf:
        writer = csv.DictWriter(outf, fieldnames = OUT_COLS, delimiter = "\t")
        writer.writeheader()
        for row in out_rows.values():
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", help = "A descriptor for the input file, used for naming output", required = True)
    parser.add_argument("--dist_report", help = "The output TSV from skani dist", required = True)
    parser.add_argument("--ani_thres", type = float, help = "Threshold for ANI", required = True)
    parser.add_argument("--align_ref_thres", type = float, help = "Threshold for % ref covered by query", required = True)
    args = parser.parse_args()

    select_refs(args.sample, args.dist_report, args.ani_thres, args.align_ref_thres)

