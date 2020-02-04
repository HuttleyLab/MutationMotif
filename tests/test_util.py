from unittest import TestCase, main

from numpy import array
from numpy.testing import assert_array_equal

from cogent3 import DNA, load_aligned_seqs
from mutation_motif.util import array_to_str, just_nucs, seqs_to_array


class TestJustNucs(TestCase):
    aln = load_aligned_seqs("data/just_nuc.fasta", array_align=True, moltype=DNA)

    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,
        seq1, seq6 to seq8 are expected to be eliminated by the code
        """
        nucs = just_nucs(self.aln.array_seqs)

        assert_array_equal(
            nucs,
            [
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0],
            ],
        )


class TestAlignSnpAnnotation(TestCase):
    d_aln = load_aligned_seqs(
        "data/load_seqs_to_array.fasta", array_align=True, moltype=DNA
    )

    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,

        expect to convert seq0, seq2 to seq5 into numpy array
        """
        data = seqs_to_array(self.d_aln)
        expect = array(
            [
                (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
            ]
        )

        assert_array_equal(data, expect)


class TestAlignSnpAnnotation2(TestCase):
    seq_array = array([(2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 0)])

    def test_array_to_str(self):
        """convert numpy array back to DNA sequence"""
        dna_str = array_to_str(self.seq_array)
        expect = ["ATCAACATATAAAAAGGAAAT"]

        self.assertEqual(dna_str, expect)


if __name__ == "__main__":
    main()
