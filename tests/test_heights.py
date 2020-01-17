from unittest import TestCase, main

from numpy.testing import assert_array_equal

from mutation_motif.height import get_mi_char_heights, get_re_char_heights
from numpy import array


class TestCalcCharHeights(TestCase):
    def test_mi_char_heigths(self):
        """calc correct char heights for MI variant"""
        probs = array(
            [
                [0.21428571, 0.07142857, 0.28571429],
                [0.5, 0.78571429, 0.28571429],
                [0.14285714, 0.07142857, 0.21428571],
                [0.14285714, 0.07142857, 0.21428571],
            ]
        )

        mi = array([0.22167, 0.91078, 0.01476])

        expect_char_mi_height = array(
            [
                [0.0475, 0, 0.00422],
                [0.11084, 0, 0.00422],
                [0.03167, 0, 0.00316],
                [0.03167, 0, 0.00316],
            ]
        )

        got_mi_height = get_mi_char_heights(probs, mi).round(decimals=5)

        assert_array_equal(got_mi_height, expect_char_mi_height)

    def test_re_char_heights(self):
        """char heights for RE terms should be correct"""
        rets = array(
            [
                [0, 0, 0.12535],
                [0.15782, 0, 0.67927],
                [0, 0, 0.28571],
                [-0.1, 0, -0.21429],
            ]
        )

        expect_char_ret_height = array(
            [
                [0, 0, 0.08417],
                [0.03539, 0, 0.45612],
                [0, 0, 0.19185],
                [-0.02243, 0, -0.14389],
            ]
        )
        got_ret_height = get_re_char_heights(rets).round(decimals=5)
        assert_array_equal(got_ret_height, expect_char_ret_height)


if __name__ == "__main__":
    main()
