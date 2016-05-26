#!/usr/bin/env python2.7
from cogent import LoadTable
from cogent.util.unit_test import TestCase, main
from mutation_motif.complement import _reverse_complement

class TestEntropy(TestCase):
    def test_reverse_complement(self):
        data = [[1670, 'T', 'T', 'T', 'T', 'M', 'AtoC'],
                [557, 'G', 'A', 'A', 'C', 'M', 'AtoC'],
                [1479, 'A', 'G', 'A', 'A', 'M', 'AtoC'],
                [925, 'G', 'A', 'A', 'G', 'M', 'AtoC'],
                [1919, 'A', 'C', 'A', 'A', 'M', 'AtoC'],
                [442, 'G', 'A', 'C', 'A', 'M', 'AtoC']]
        header = ['count', 'pos0', 'pos1', 'pos2', 'pos3', 'mut', 'direction']
        table = LoadTable(header=header, rows=data)
        ex =  [[1670, 'A', 'A', 'A', 'A', 'M', 'AtoC'],
                [557, 'G', 'T', 'T', 'C', 'M', 'AtoC'],
                [1479, 'T', 'T', 'C', 'T', 'M', 'AtoC'],
                [925, 'C', 'T', 'T', 'C', 'M', 'AtoC'],
                [1919, 'T', 'T', 'G', 'T', 'M', 'AtoC'],
                [442, 'T', 'G', 'T', 'C', 'M', 'AtoC']]
        got = _reverse_complement(table)
        raw_got = got.getRawData()
        
        self.assertEqual(raw_got, ex)

if __name__ == '__main__':
    main()

