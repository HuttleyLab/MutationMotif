#!/usr/bin/env python2.7
from cogent import LoadTable, DNA
from cogent.util.unit_test import TestCase, main
from mutation_motif.complement import _reverse_complement,\
                    make_strand_symmetric_table

class TestEntropy(TestCase):
    data = [[1670, 'T', 'T', 'T', 'T', 'M', 'AtoC'],
            [557, 'G', 'A', 'A', 'C', 'M', 'AtoC'],
            [1479, 'A', 'G', 'A', 'A', 'M', 'AtoC'],
            [925, 'G', 'A', 'A', 'G', 'M', 'AtoC'],
            [1919, 'A', 'C', 'A', 'A', 'M', 'AtoC'],
            [442, 'G', 'A', 'C', 'A', 'M', 'AtoC']]
    header = ['count', 'pos0', 'pos1', 'pos2', 'pos3', 'mut', 'direction']
    
    def test_reverse_complement(self):
        table = LoadTable(header=self.header, rows=self.data)
        ex =  [[1670, 'A', 'A', 'A', 'A', 'M', 'AtoC'],
                [557, 'G', 'T', 'T', 'C', 'M', 'AtoC'],
                [1479, 'T', 'T', 'C', 'T', 'M', 'AtoC'],
                [925, 'C', 'T', 'T', 'C', 'M', 'AtoC'],
                [1919, 'T', 'T', 'G', 'T', 'M', 'AtoC'],
                [442, 'T', 'G', 'T', 'C', 'M', 'AtoC']]
        got = _reverse_complement(table)
        raw_got = got.getRawData()
        
        self.assertEqual(raw_got, ex)
    
    def test_strandsym_table(self):
        """makes strand symmetric table"""
        data = [[1, 'T', 'T', 'T', 'T', 'M', 'TtoG'],
                [1, 'G', 'A', 'A', 'C', 'M', 'TtoG'],
                [1, 'A', 'G', 'A', 'A', 'M', 'TtoG'],
                [1, 'G', 'A', 'A', 'G', 'M', 'TtoG'],
                [1, 'A', 'C', 'A', 'A', 'M', 'TtoG'],
                [1, 'G', 'A', 'C', 'A', 'M', 'TtoG']]
        exp = []
        for row in self.data:
            n = row[:]
            n.append('+')
            exp.append(n)
        for row in data:
            seq = map(DNA.complement, row[1:-2])
            seq.reverse()
            n = [row[0]] + seq + ['M', 'AtoC']
            n.append('-')
            exp.append(n)
        
        table = LoadTable(header=self.header, rows=self.data+data)
        r = make_strand_symmetric_table(table)
        self.assertEqual(r.getRawData(), exp)

if __name__ == '__main__':
    main()

