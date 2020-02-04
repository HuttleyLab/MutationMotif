from unittest import TestCase, main

from numpy import array
from numpy.testing import assert_array_equal, assert_allclose

from mutation_motif.height import get_mi_char_heights, get_re_char_heights


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

        expect_char_mi_height = probs * mi

        got_mi_height = get_mi_char_heights(probs, mi)

        assert_allclose(got_mi_height, expect_char_mi_height)

    def test_re_char_heights(self):
        """char heights for RE terms should be correct"""
        rets = array(
            [
                [0,      0,  0.12535],
                [0.15782, 0, 0.67927],
                [0,      0,  0.28571],
                [-0.1,   0,  -0.21429],
            ]
        )

        expect_char_ret_height = array(
            [
                [0,       0, 0.08417],
                [0.035394, 0, 0.45612],
                [0,       0, 0.19185],
                [-0.022426, 0, -0.14389],
            ]
        )
        got_ret_height = get_re_char_heights(rets)
        assert_allclose(got_ret_height, expect_char_ret_height, rtol=1e-4)


if __name__ == "__main__":
    main()
