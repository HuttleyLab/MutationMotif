from unittest import TestCase, main

from cogent3 import DNA
from cogent3.core.alignment import ArrayAlignment

from mutation_motif.motif_count import (
    get_combined_counts,
    get_count_table,
    profile_to_seq_counts,
    reduced_multiple_positions,
    reduced_one_position,
)


obs = """
TATGT
TATTT
CTTTA
TATCA
TATCA
CTTCT
""".splitlines()

obs = [("s%s" % i, s) for i, s in enumerate(obs) if s]

ctl = """
TCTTC
TCTTC
GGTAA
AGTGT
GCTGG
GGTAG
""".splitlines()

ctl = [("s%s" % i, s) for i, s in enumerate(ctl) if s]


def _get_seq_array(data):
    """returns [(n, seq), ...] as DenseArray"""
    return ArrayAlignment(data=data, moltype=DNA).array_seqs


class TestMotifCount(TestCase):
    obs_count = profile_to_seq_counts(_get_seq_array(obs), 2)
    ctl_count = profile_to_seq_counts(_get_seq_array(ctl), 2)
    table = get_count_table(obs_count, ctl_count, 4)

    def test_profile_to_counts(self):
        """correctly determine counts of motif flanking central base"""
        seqs = _get_seq_array(obs)
        result = profile_to_seq_counts(seqs, 2)
        exp = dict(TAGT=1, TATT=1, CTTA=1, TACA=2, CTCT=1)
        self.assertEqual(result, exp)

    def test_counts_table(self):
        """construct counts table"""
        obs_count = profile_to_seq_counts(_get_seq_array(obs), 2)
        ctl_count = profile_to_seq_counts(_get_seq_array(ctl), 2)

        r = get_count_table(obs_count, ctl_count, 4)
        self.assertEqual(r.distinct_values("mut"), set("MR"))
        # because the motifs are unique in ctl and obs
        # total number should double
        self.assertEqual(r.shape[0], 2 * (len(obs_count) + len(ctl_count)))

    def test_reduced(self):
        """reduced across positions should produce correct counts"""
        r = reduced_multiple_positions(self.table, "pos0")
        exp = {
            "R": {("A",): 1, ("C",): 0, ("G",): 3, ("T",): 2},
            "M": {("A",): 0, ("C",): 2, ("G",): 0, ("T",): 4},
        }
        self.assertEqual(r, exp)

        r = reduced_multiple_positions(self.table, "pos2")
        exp = {
            "R": {("A",): 2, ("C",): 0, ("G",): 2, ("T",): 2},
            "M": {("A",): 0, ("C",): 3, ("G",): 1, ("T",): 2},
        }
        self.assertEqual(r, exp)

        r = reduced_multiple_positions(self.table, "pos0", "pos2")
        exp = {
            "R": {
                ("A", "G"): 1,
                ("C", "C"): 0,
                ("C", "T"): 0,
                ("G", "G"): 1,
                ("G", "A"): 2,
                ("T", "C"): 0,
                ("T", "G"): 0,
                ("T", "T"): 2,
            },
            "M": {
                ("A", "G"): 0,
                ("C", "C"): 1,
                ("C", "T"): 1,
                ("G", "G"): 0,
                ("G", "A"): 0,
                ("T", "C"): 2,
                ("T", "G"): 1,
                ("T", "T"): 1,
            },
        }
        self.assertEqual(r, exp)

    def test_reduced_one(self):
        """should give same results as more general"""
        r = reduced_multiple_positions(self.table, "pos0")
        got = reduced_one_position(self.table, "pos0")
        exp = {}
        for m, counts in list(r.items()):
            counts = dict((k[0], c) for k, c in list(counts.items()))
            exp[m] = counts
        self.assertEqual(got, exp)

        r = reduced_multiple_positions(self.table, "pos2")
        exp = {}
        for m, counts in list(r.items()):
            counts = dict((k[0], c) for k, c in list(counts.items()))
            exp[m] = counts

        got = reduced_one_position(self.table, "pos2")
        self.assertEqual(got, exp)

    def test_combined_counts(self):
        """combining counts completes missing states"""
        combined = get_combined_counts(self.table, "pos0")
        self.assertEqual(combined.shape[0], 8)
        combined = get_combined_counts(self.table, ["pos0", "pos1"])
        self.assertEqual(combined.shape[0], 32)

        # the following states are not present in either group for pos0/1
        # and so should be missing
        missing = [
            ["A", "A"],
            ["A", "C"],
            ["A", "T"],
            ["C", "C"],
            ["C", "G"],
            ["C", "A"],
            ["G", "A"],
            ["G", "T"],
            ["T", "G"],
            ["T", "T"],
        ]
        for b1, b2 in missing:
            sub = combined.filtered("base1=='%s' and base2=='%s'" % (b1, b2))
            self.assertEqual(sub.summed("count"), 0)


if __name__ == "__main__":
    main()
