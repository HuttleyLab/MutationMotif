#!/usr/bin/env python2.7
import sys, os
from cogent import LoadSeqs, DNA
from cogent.util.unit_test import TestCase, main

from mutation_motif.height import get_re_char_heights, get_mi_char_heights

class TestCalcCharHeights(TestCase):
    
    def test_re_char_heights(self):
        """char heights for RE terms should be correct
        """
        # make an example with length 3 with specified values for the RETs
        # e.g. [[-0.4, -0.2, -0.1, 0.4], [0,0,0,0],# SNP loc
        #        [-0.1, 0.1, 0.2, 0.3]]
        # from these you should know the bar height and the proportions
        raise NotImplementedError
    
    def test_mi_char_heigths(self):
        """calc correct char heights for MI variant
        """
        raise NotImplementedError
    
if __name__ == '__main__':
    main()

        