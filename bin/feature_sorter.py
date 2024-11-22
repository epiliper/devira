import itertools


def pairwise(iterable):
    """from itertools recipes
    s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


class FeatureSorter(object):
    """This class helps sort genomic features. It's not terribly optimized
    for speed or anything. Slightly inspired by calhoun's MultiSequenceRangeMap.
    """

    def __init__(self, collection=None):
        self.seqids = []
        self.seq_to_features = {}
        self.seq_to_breakpoints = {}
        self.dirty = False
        if collection is not None:
            for args in collection:
                self.add(*args)

    def add(self, c, start, stop, strand="+", other=None):
        """Add a "feature", which is a chrom,start,stop,strand tuple (with
        optional other info attached)
        """
        if c not in self.seq_to_features:
            self.seqids.append(c)
            self.seq_to_features[c] = []
            self.seq_to_breakpoints[c] = set()
            # self.seq_to_breakpoints[c].add(1) # do we want the empty interval in front?
        self.seq_to_features[c].append((int(start), int(stop), strand, other))
        self.seq_to_breakpoints[c].add(start)
        self.seq_to_breakpoints[c].add(stop + 1)
        self.dirty = True

    def _cleanup(self):
        if self.dirty:
            self.dirty = False
            for c in self.seqids:
                self.seq_to_features[c].sort()

    def get_seqids(self):
        return tuple(self.seqids)

    def get_features(self, c=None, left=0, right=float("inf")):
        """Get all features on all chromosomes in sorted order. Chromosomes
        are emitted in order of first appearance (via add). Features on
        each chromosome are emitted in start, then stop order.  If
        boundaries are specified, we restrict to features that contain
        the specified interval.
        """
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for start, stop, strand, other in self.seq_to_features[c]:
                if stop >= left and start <= right:
                    yield (c, start, stop, strand, other)

    def get_intervals(self, c=None):
        """Get all intervals on the reference where the overlapping feature
        set remains the same. Output will be sorted, adjacent intervals
        and will describe how many and which features overlap it.
        """
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for left, right in pairwise(sorted(self.seq_to_breakpoints[c])):
                right = right - 1
                features = list(self.get_features(c, left, right))
                yield (c, left, right, len(features), features)
