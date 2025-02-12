#!/usr/bin/env python3

# note: this code is a simpler version of
# https://github.com/broadinstitute/viral-assemble/blob/master/assemble/mummer.py
# prone to change
# so far just works with single-reference (but multi-segment) fastas as inputs
# also, no support for looking at alternate contigs/sequences

from feature_sorter import FeatureSorter
from aligns_reader import AlignsReader
import random
import subprocess
import logging
import Bio.SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tile_file")
parser.add_argument("ref_fasta")
parser.add_argument("delta_file")
parser.add_argument("--ambig_max_aligns", default=2, required=False)
parser.add_argument("--ambig_max_lens", default=1, required=False)
parser.add_argument("--ambig_max_frac", default=0.01, required=False)
parser.add_argument("--out_scaffold_name")
args = parser.parse_args()


log = logging.getLogger(__name__)


def histogram(items):
    """I count the number of times I see stuff and return a dict of counts."""
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
    return out


def fastaMaker(seqs, linewidth=60):
    assert linewidth > 0

    for idVal, seq in seqs:
        yield ">{}\n".format(idVal)

        while len(seq) > linewidth:
            line = seq[:linewidth]
            seq = seq[linewidth:]
            yield "{}\n".format(line)

        if seq:
            yield seq + "\n"


fs = FeatureSorter()

def contig_chooser(alt_seqs, ref_len, coords_debug=""):
    """Our little heuristic to choose an alternative sequence from a pile
    of alignments to a reference. Takes a list of strings (one string per
    contig). This method will choose a single sequence from the input:
        1. if there are no alt_seqs, emit a stretch of Ns, same length as ref
        2. if there is only one alt_seq, emit that one
        3. if there are many alt_seqs, emit the most popular (if there is one)
        4. otherwise, if there is a most popular sequence length (including
            the ref_len as one vote), filter the alt_seqs to those of the
            most popular length and emit the most popular sequence (if there
            is one), otherwise, choose randomly amongst the remainder
        5. otherwise, choose randomly amongst the same-as-ref length sequences
        6. or just choose randomly if there are no same-as-ref length sequences
    The output will be a list of unique strings, where the first string
    is the "chosen" sequence, and the remaining strings are the alternative
    sequences (in no particular order, but guaranteed to be unique).
    """
    if not alt_seqs:
        # no contigs align here, emit Ns of appropriate length
        new_seq = "N" * ref_len
        other_seqs = []
    elif len(alt_seqs) == 1:
        # only one contig aligns here
        new_seq = alt_seqs[0]
        other_seqs = []
    else:
        # multiple contigs align here
        ranks = list(
            sorted(histogram(alt_seqs).items(), key=lambda x: x[1], reverse=True)
        )
        other_seqs = list(s for s, n in ranks)
        if len(ranks) == 1:
            # only one unique sequence exists
            new_seq = ranks[0][0]
        elif ranks[0][1] > ranks[1][1]:
            # clear winner: a most popular sequence exists
            new_seq = ranks[0][0]
        else:
            # multiple possible replacement sequences
            len_ranks = list(
                sorted(
                    histogram(
                        [len(s) for s in alt_seqs]
                        + [ref_len]  # let the reference have one vote
                    ).items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
            )
            if len(len_ranks) == 1 or len_ranks[0][1] > len_ranks[1][1]:
                # a most popular replacement length exists
                # remove all alt_seqs that don't use that length
                alt_seqs = list(s for s in alt_seqs if len(s) == len_ranks[0][0])
                assert alt_seqs
                ranks = list(
                    sorted(
                        histogram(alt_seqs).items(),
                        key=lambda x: x[1],
                        reverse=True,
                    )
                )
                if len(ranks) == 1 or ranks[0][1] > ranks[1][1]:
                    # clear winner amongst remaining sequences of most popular length
                    new_seq = ranks[0][0]
                else:
                    # more complicated scenario. choose randomly.
                    # perhaps in future, vote based on aligned read count?
                    if len(alt_seqs) > 1:
                        log.warning(
                            "choosing random contig from %d choices of most popular length in %s"
                            % (len(alt_seqs), coords_debug)
                        )
                    new_seq = random.choice(alt_seqs)
            else:
                # no clear winner on replacement length
                alt_ref_len_seqs = list(s for s in alt_seqs if len(s) == ref_len)
                if alt_ref_len_seqs:
                    # choose randomly among same-as-ref-length sequences
                    alt_seqs = alt_ref_len_seqs
                    if len(alt_seqs) > 1:
                        log.warning(
                            "choosing random contig from %d choices of reference length in %s"
                            % (len(alt_seqs), coords_debug)
                        )
                    new_seq = random.choice(alt_seqs)
                else:
                    # no clear winner and all replacement lengths are different from reference length
                    # just choose randomly
                    if len(alt_seqs) > 1:
                        log.warning(
                            "choosing random contig from %d choices in %s"
                            % (len(alt_seqs), coords_debug)
                        )
                    new_seq = random.choice(alt_seqs)
        other_seqs = list(s for s in other_seqs if s != new_seq)
    return [new_seq] + other_seqs


def get_best_from_tiling(
    tile_file,
    ref_file,
    delta_file,
    ambig_max_lens,
    ambig_max_aligns,
    ambig_max_frac,
    out_scaffolds_fasta,
):
    with open(tile_file, "rt") as tilef:
        for line in tilef:
            row = line.rstrip("\n\r").split("\t")
            c = row[11]
            start, stop = (int(row[0]), int(row[1]))
            alt_seq = (row[12], int(row[2]), int(row[3]))
            if stop < start:
                raise ValueError()
            if alt_seq[2] < alt_seq[1]:
                s = "-"
            else:
                s = "+"

            fs.add(c, start, stop, strand=s, other=alt_seq)

    alnReaders = {}

    aln_file = f"{random.randint(100000, 999999)}.aligns"

    for c, start, stop, strand, other in fs.get_features():
        chr_pair = (c, other[0])
        if chr_pair not in alnReaders:
            toolCmd = [
                # "mummer",
                "show-aligns",
                "-r",
                delta_file,
                chr_pair[0],
                chr_pair[1],
            ]
            with open(aln_file, "wt") as outf:
                subprocess.check_call(toolCmd, stdout=outf)
            alnReaders[chr_pair] = AlignsReader(aln_file, ref_file)

    with open(out_scaffolds_fasta, "wt") as outf:
        for c in [seqObj.id for seqObj in Bio.SeqIO.parse(ref_file, "fasta")]:
            if c not in fs.get_seqids():
                # old: continue
                log.FATAL("ref for contig not found")

            def n_diff_vals(*vals):
                return len(set(vals))

            def n_diff_lens(seqs):
                return n_diff_vals(*map(len, seqs))

            def frac_unambig(seqs):
                """Given a list of seqs of the same length, return the fraction of positions on which they all agree"""
                n_tot = len(seqs[0])
                n_unambig = list(map(n_diff_vals, *seqs)).count(1)
                return float(n_unambig) / float(n_tot or 1.0)

            # construct scaffolded sequence for this chromosome
            seq = []
            for _, left, right, n_features, features in fs.get_intervals(c):
                # get all proposed sequences for this specific region
                alt_seqs = []
                for consider_ambig_aligns in (False, True):
                    for f in features:
                        alt_seqs_f = alnReaders[
                            (c, f[-1][0])
                        ].retrieve_alts_by_ref(
                            left, right, aln_start=f[1], aln_stop=f[2]
                        )
                        if len(alt_seqs_f) == 1:
                            alt_seqs.append(alt_seqs_f[0])
                        elif consider_ambig_aligns:
                            if (
                                len(alt_seqs_f) <= ambig_max_aligns
                                and n_diff_lens(alt_seqs_f) <= ambig_max_lens
                                and frac_unambig(alt_seqs_f)
                                > (1.0 - ambig_max_frac)
                            ):
                                alt_seqs.append(alt_seqs_f[0])
                                log.info(
                                    "using ambiguous alignment to ref seq {} at [{},{}]".format(
                                        c, f[1], f[2]
                                    )
                                )
                            else:
                                log.warning(
                                    "dropping ambiguous alignment to ref seq {} at [{},{}]".format(
                                        c, f[1], f[2]
                                    )
                                )
                    if alt_seqs:
                        # if have a non-unambiguous alignment, don't consider ambiguous ones
                        break

                # pick the "right" one and glue together into a chromosome
                ranked_unique_seqs = contig_chooser(
                    alt_seqs, right - left + 1, "%s:%d-%d" % (c, left, right)
                )
                seq.append(ranked_unique_seqs[0])

        for line in fastaMaker(
            [(str(c) + "_contigs_ordered_and_oriented", "".join(seq))]
        ):
            outf.write(line)

if __name__ == "__main__":
    get_best_from_tiling(
        args.tile_file,
        args.ref_fasta,
        args.delta_file,
        args.ambig_max_lens,
        args.ambig_max_aligns,
        args.ambig_max_frac,
        args.out_scaffold_name,
    )
