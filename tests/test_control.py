from unittest import TestCase, main

from numpy import array
from numpy.testing import assert_array_equal

from cogent3 import DNA, make_aligned_seqs
from mutation_motif.profile import (
    MakeCircleRange,
    chosen_base_indices,
    filter_seqs_by_chosen_base,
    get_control,
    get_random_indices,
    get_zero_counts,
)


class TestChooseBases(TestCase):
    data_even = array(
        [
            (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3),
        ]
    )

    data_odd = array(
        [
            (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
        ]
    )

    chosen_pair = "AC"
    chosen_base = "A"  # 2

    step1 = 1
    step2 = 3

    def test_even_length_1(self):
        """raise the AssertionError when seqs not 2n + 1 long or chosen base is not a single base"""
        self.assertRaises(
            AssertionError,
            chosen_base_indices,
            self.data_even,
            self.chosen_base,
            self.step1,
        )
        self.assertRaises(
            AssertionError,
            chosen_base_indices,
            self.data_even,
            self.chosen_base,
            self.step2,
        )
        self.assertRaises(
            AssertionError,
            chosen_base_indices,
            self.data_odd,
            self.chosen_pair,
            self.step1,
        )
        self.assertRaises(
            AssertionError,
            chosen_base_indices,
            self.data_odd,
            self.chosen_pair,
            self.step2,
        )

    def test_chosen_base_indices1(self):
        """when setting the picking interval as 1,
        this piece of code should produce a list of
        potential pseudo-SNP loc-indicies
        """
        indicies = chosen_base_indices(self.data_odd, self.chosen_base, self.step1)
        expect = [
            array(
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            ),  # every position valid
            array([]),  # no valid positions
            array([0, 3, 6, 7, 8, 11, 12, 13, 14, 17, 18]),
        ]  # some valid positions

        for i in range(len(indicies)):
            assert_array_equal(indicies[i], expect[i])

    def test_chosen_base_indices3(self):
        """when setting the picking interval as 3,
        this piece of code should produce a list of
        potential pseudo-SNP loc-indicies
        """
        indicies = chosen_base_indices(self.data_odd, self.chosen_base, self.step2)
        expect = [
            array([1, 4, 7, 13, 16, 19]),  # every position valid
            array([]),  # no valid positions
            array([7, 13]),
        ]  # some valid positions

        for i in range(len(indicies)):
            assert_array_equal(indicies[i], expect[i])


class TestAlignFiltering(TestCase):
    data = array(
        [
            (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
        ]
    )

    sampled_indices = [array([1, 4, 7, 13, 16, 19]), array([]), array([7, 13])]

    min_chosen_bases1 = 1
    min_chosen_bases2 = 4

    def test_filter_seqs_1(self):
        """only seqs with >= 1 potential pseudo-SNP base should be returned
        """
        test_data, sample_indices = filter_seqs_by_chosen_base(
            self.data, self.sampled_indices, self.min_chosen_bases1
        )

        assert_array_equal(
            test_data,
            array(
                [
                    (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                    (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
                ]
            ),
        )

        assert_array_equal(sample_indices[0], array([1, 4, 7, 13, 16, 19]))
        assert_array_equal(sample_indices[1], array([7, 13]))

    def test_filter_seqs_4(self):
        """only seqs with >= 4 potential pseudo-SNP base should be returned
        """
        test_data, sample_indices = filter_seqs_by_chosen_base(
            self.data, self.sampled_indices, self.min_chosen_bases2
        )

        assert_array_equal(
            test_data,
            array([(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)]),
        )

        assert_array_equal(sample_indices[0], array([1, 4, 7, 13, 16, 19]))


class TestRandomIndices(TestCase):
    seq_array = array([(2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 0)])

    sample_indices = [array([4, 13, 19])]

    slice_side = 7

    def test_get_random_indices(self):
        """the testing seq array contains 3 potential pseudo-snp,
        the produced list should be one of the 3 possible potential
        pseudo-snps and its coordinate flanking seq"""
        expected = [
            [[18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]],
            [[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]],
            [[12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5]],
        ]
        circle_range = MakeCircleRange(self.seq_array.shape[1], self.slice_side)
        for i in range(4):
            sampled_indices = get_random_indices(self.sample_indices, circle_range)
            self.assertTrue(sampled_indices in expected)


class TestAlignSnpAnnotation(TestCase):
    seqs = [("seq_0", "ATCAACATATAAAAAGGAAAT")]

    d_aln = make_aligned_seqs(data=seqs, array_align=True, moltype=DNA).array_seqs
    step = 3
    slice_side = 7
    direction = "AtoC"
    seed = "1234"
    chosen_base = direction[0]

    def test_get_control(self):
        """return a valid control
        """
        expected = [
            [[2, 2, 0, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 2, 2]],
            [[2, 0, 2, 0, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 0]],
            [[2, 2, 2, 3, 3, 2, 2, 2, 0, 2, 0, 1, 2, 2, 1]],
        ]

        for i in range(5):
            control = get_control(
                self.d_aln, self.chosen_base, self.step, self.slice_side, seed=self.seed
            )
            self.assertTrue(control.tolist() in expected)


class TestIndices(TestCase):
    circle_size_1 = 20
    circle_size_2 = 21
    slice_side = 7

    def test_even_circle_size(self):
        """raise the AssertionError when circle size is not 2n + 1 long"""
        self.assertRaises(
            AssertionError, MakeCircleRange, self.circle_size_1, self.slice_side
        )

    def test_MakeCircleRange(self):
        """this test should prduce a full set of possible circular slices,
        slices are assumed to be 2*7+1=15 in length
        """
        circle_range = MakeCircleRange(self.circle_size_2, self.slice_side)
        full_list = []

        for i in range(self.circle_size_2):
            full_list.append(circle_range(i))

        expect = [
            [14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7],
            [15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8],
            [16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            [17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            [20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
            [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
            [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0],
            [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1],
            [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2],
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3],
            [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4],
            [12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5],
            [13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6],
        ]

        self.assertEqual(full_list, expect)


class TestProfile(TestCase):
    def test_zero_counts(self):
        """zero profile constructed correctly"""
        c = get_zero_counts(5, int, pseudo_count=0)
        self.assertEqual(c.shape, (4, 5))
        self.assertEqual(c.tolist(), [[0, 0, 0, 0, 0]] * 4)
        c = get_zero_counts(5, int, pseudo_count=1)
        self.assertEqual(c.tolist(), [[1, 1, 1, 1, 1]] * 4)


if __name__ == "__main__":
    main()
