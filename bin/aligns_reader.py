import logging
import Bio


log = logging.getLogger(__name__)


class AlignsReader(object):
    """This class assists in the parsing and reading of show-aligns output."""

    def __init__(self, aligns_file, ref_fasta=None):
        self.aligns_file = aligns_file
        self.ref_fasta = ref_fasta
        self.reference_seq = None
        self.alignments = []
        self.seq_ids = []
        self._load_align()
        self._load_fastas()

    def _load_align(self):
        with open(self.aligns_file, "rt") as inf:
            # read ref_fasta, query_fasta from header
            header = inf.readline().strip().split()
            assert len(header) == 2
            if self.ref_fasta is None:
                self.ref_fasta = header[0]

            # iterate row by row
            mode = "start"
            coords = None
            seqs = None
            for line in inf:
                line = line.rstrip()
                if not line:
                    pass  # empty lines
                elif not line.strip("="):
                    pass  # header and footer of file
                elif line[0] in (" ", "\t"):
                    pass  # describes mismatches in alignment
                elif line.startswith("--"):
                    if line.startswith("-- Alignments between "):
                        assert mode == "start"
                        mode = "between"
                        self.seq_ids = line[22:].split(" and ")
                        assert len(self.seq_ids) == 2
                    elif line.startswith("-- BEGIN alignment [ "):
                        assert mode == "between"
                        mode = "align"
                        coords = list(x.split() for x in line[21:-2].split(" | "))
                        assert len(coords) == 2 and coords[0][2] == coords[1][2] == "-"
                        assert coords[0][0] == "+1", "error with line: %s" % line
                        seqs = [[], []]
                        align_lines = 0
                    elif line.startswith("--   END alignment [ "):
                        assert mode == "align"
                        mode = "between"
                        new_coords = list(x.split() for x in line[21:-2].split(" | "))
                        assert coords == new_coords, "error: %s != %s" % (
                            new_coords,
                            coords,
                        )
                        assert len(seqs[0]) == len(seqs[1])
                        seqs = list("".join(x) for x in seqs)
                        assert len(seqs[0]) == len(seqs[1])
                        self.alignments.append(
                            [
                                coords[0][0],
                                int(coords[0][1]),
                                int(coords[0][3]),
                                coords[1][0],
                                int(coords[1][1]),
                                int(coords[1][3]),
                                seqs[0],
                                seqs[1],
                            ]
                        )
                        coords = None
                        seqs = None
                    else:
                        raise AssertionError("file format: line '%s'" % line)
                else:
                    # read in one line of an alignment
                    assert mode == "align", (
                        "file format: line '%s' before alignment begins" % line
                    )
                    seq_str = line.split()[1].upper().replace(".", "-")
                    seqs[align_lines % 2].append(seq_str)
                    align_lines += 1

    def _load_fastas(self):
        assert self.ref_fasta and self.seq_ids
        self.reference_seq = Bio.SeqIO.index(self.ref_fasta, "fasta")[self.seq_ids[0]]

    def get_alignments(self):
        for a in self.alignments:
            yield a

    def get_intervals(self):
        prev = None
        for a in self.alignments:
            cur = a[1:3]
            assert cur[1] >= cur[0]
            if prev is None:
                if cur[0] > 1:
                    # emit leading reference sequence before first alignment
                    yield self._dummy_row(1, cur[0] - 1, "-")
            else:
                assert cur[0] > prev[1], "overlaps not allowed"
                if cur[0] > prev[1] + 1:
                    # emit gap between alignments
                    yield self._dummy_row(prev[1] + 1, cur[0] - 1, "N")
            # emit actual alignment
            yield a
            prev = cur
        if prev and prev[1] < len(self.reference_seq):
            # emit trailing reference sequence after last alignment
            yield self._dummy_row(prev[1] + 1, len(self.reference_seq), "-")

    def _dummy_row(self, start, stop, filler="N"):
        return [
            "+1",
            start,
            stop,
            "+1",
            start,
            stop,
            self.get_ref_seq(start, stop),
            filler * (stop - start + 1),
        ]

    def get_ref_seq(self, start, stop):
        """Retrieve a sub-sequence from the reference (1st) sequence in the
        alignment using coordinates relative to the reference sequence.
        No gaps will be emitted.
        """
        return str(self.reference_seq.seq[start - 1 : stop])

    def retrieve_alts_by_ref(self, start, stop, aln_start=None, aln_stop=None):
        """Retrieve sub-sequence(s) from the alternate (2nd) sequence in the
        alignment using coordinates relative to the reference sequence.
        No gaps will be emitted.
        Required: start-stop interval must be wholly contained within
        an alignment.
        """

        # grab the one alignment that contains this window
        alns = list(a for a in self.alignments if a[1] <= start and a[2] >= stop)
        if aln_start is not None and aln_stop is not None:
            # if specified, restrict to a specific alignment that comes from show-tiling
            # (sometimes show-aligns is more promiscuous than show-tiling)
            new_alns = []
            for a in alns:
                if a[1] > aln_start or a[2] < aln_stop:
                    log.debug(
                        "dropping undesired alignment: %s(%s):%s-%s to %s(%s):%s-%s (%s:%s-%s requested)",
                        self.seq_ids[0],
                        a[0],
                        a[1],
                        a[2],
                        self.seq_ids[1],
                        a[3],
                        a[4],
                        a[5],
                        self.seq_ids[0],
                        aln_start,
                        aln_stop,
                    )
                else:
                    new_alns.append(a)
            alns = new_alns
        if len(alns) != 1:
            log.warning(
                "invalid %s:%d-%d -> %s specified, %d alignments found that contain it",
                self.seq_ids[0],
                start,
                stop,
                self.seq_ids[1],
                len(alns),
            )
            for aln in alns:
                log.debug("alignment: %s", str(aln[:6]))

        return [self._aln_to_alt_seq(aln, start, stop) for aln in alns]

    def _aln_to_alt_seq(self, aln, start, stop):
        """Given an alignment of a contig to ref, return the contig sequence aligned to a given stretch of ref"""
        ref_l, ref_r, ref_seq, alt_seq = (aln[1], aln[2], aln[-2], aln[-1])

        # convert desired start/stop relative to this reference window
        #  such that 0 <= start <= stop <= ref_r-ref_l+1
        aln_start = start - ref_l
        aln_stop = stop - ref_l

        # travel down alignment until we've reached the left edge
        #  (because of gaps, you must check each position one by one)
        #  end loop when ref_seq[:i_left] contains {aln_start} bases
        n_ref_bases = 0
        i_left = 0
        while n_ref_bases < aln_start:
            if ref_seq[i_left] != "-":
                n_ref_bases += 1
            i_left += 1

        # travel down alignment until we've reached the right edge
        #  (because of gaps, you must check each position one by one)
        #  end loop when ref_seq[:i_right] contains {aln_stop} bases
        i_right = i_left
        while n_ref_bases < aln_stop:
            if ref_seq[i_right] != "-":
                n_ref_bases += 1
            i_right += 1
        # consume and include any trailing gaps
        while i_right < len(ref_seq) and ref_seq[i_right] == "-":
            i_right += 1

        # grab the alternate sequence and strip gaps
        return alt_seq[i_left : i_right + 1].replace("-", "")
