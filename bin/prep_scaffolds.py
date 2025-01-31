#!/usr/bin/env python3

from aligns_reader import AlignsReader
import Bio.AlignIO
import Bio.SeqIO
import argparse
import subprocess
import random

parser = argparse.ArgumentParser()
parser.add_argument("delta_file")
parser.add_argument("query_fasta")
parser.add_argument("ref_fasta")
parser.add_argument("output")
args = parser.parse_args()


def scaffolds_with_reference(delta_file, query_fasta, ref_fasta):

    ref_id = Bio.SeqIO.read(ref_fasta, 'fasta').id
    query_id = Bio.SeqIO.read(query_fasta, 'fasta').id


    align_file = f"{random.randint(100000, 999999)}.aligns"

    with open(align_file, "wt") as outf:
        subprocess.check_call(["show-aligns", "-r", str(delta_file), str(ref_id), str(query_id)], stdout=outf)

    alns = AlignsReader(align_file, ref_fasta)

    seqs = [[alns.seq_ids[0], []], [alns.seq_ids[1], []]]

    for a in alns.get_intervals():
        seqs[0][1].append(a[6])
        seqs[1][1].append(a[7])

    seqs[0][1] = ''.join(seqs[0][1])
    seqs[1][1] = ''.join(seqs[1][1])

    out_fasta = f"{random.randint(100000, 999999)}.fasta"
    with open(out_fasta, "wt") as outf:

        for seq in seqs:
            outf.write(f">{seq[0]}\n")
            outf.write(f"{seq[1]}\n")

    return ref_id, query_id, out_fasta


def modify_contigs(infile, ref, query, outfile):
    aln = Bio.AlignIO.read(infile, "fasta")
    print(aln[0].name)
    print(aln[1].name)

    if len(aln) != 2:
        raise Exception("alignment does not contain exactly 2 sequences, %s found" % len(aln))
    elif aln[0].name == ref:
        ref_idx = 0
        consensus_idx = 1
    elif aln[1].name == ref:
        ref_idx = 1
        consensus_idx = 0
    else:
        raise NameError("reference name '%s' not in alignment" % ref)

    mc = ContigModifier(str(aln[ref_idx].seq), str(aln[consensus_idx].seq))
    mc.call_reference_ns()
    mc.trim_ends()
    mc.replace_5ends(55)
    mc.replace_3ends(55)
    mc.replace_end_gaps()
    mc.replace_end_gaps()

    with open(outfile, "wt") as f:
        name = aln[consensus_idx].name

        f.write(f">{name}\n")
        f.write(mc.get_stripped_consensus())

class ContigModifier(object):
    ''' Initial modifications to spades+MUMmer assembly output based on
        MUSCLE alignment to known reference genome
        author: rsealfon
    '''

    def __init__(self, ref, consensus):
        if len(ref) != len(consensus):
            raise Exception("improper alignment")
        self.ref = list(ref)
        self.consensus = list(consensus)
        self.len = len(ref)

    def get_stripped_consensus(self):
        return ''.join(self.consensus).replace('-', '')

    def call_reference_ns(self):
        for i in range(self.len):
            if self.consensus[i].upper() == "N":
                self.consensus[i] = self.ref[i]

    def call_reference_ambiguous(self):
        ''' This is not normally used by default in our pipeline '''
        for i in range(self.len):
            if self.ref[i].upper() in Bio.Data.IUPACData.ambiguous_dna_values.get(self.consensus[i].upper(), []):
                self.consensus[i] = self.ref[i]

    def trim_ends(self):
        ''' This trims down the consensus so it cannot go beyond the given reference genome '''
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.ref[i] != "-":
                    break
                else:
                    self.consensus[i] = "-"

    def replace_end_gaps(self):
        ''' This fills out the ends of the consensus with reference sequence '''
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.consensus[i] != "-":
                    break
                self.consensus[i] = self.ref[i]

    def replace_5ends(self, replace_length):
        ''' This replaces everything within <replace_length> of the ends of the
            reference genome with the reference genome.
        '''
        ct = 0
        for i in range(self.len):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i + 1):
                    self.consensus[j] = self.ref[j]
                break

    def replace_3ends(self, replace_length):
        ct = 0
        for i in reversed(range(self.len)):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i, self.len):
                    self.consensus[j] = self.ref[j]
                break

    def remove_end_ns(self):
        ''' This clips off any N's that begin or end the consensus.
            Not normally used in our pipeline
        '''
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if (self.consensus[i].upper() == "N" or self.consensus[i] == "-"):
                    self.consensus[i] = "-"
                else:
                    break

if __name__ == "__main__":
    ref_id, query_id, merged_fasta = scaffolds_with_reference(args.delta_file, args.query_fasta, args.ref_fasta)
    modify_contigs(merged_fasta, ref_id, query_id, args.output)
