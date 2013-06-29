from numpy import array

from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment
from cogent.util.unit_test import TestCase, main

from mutation_motif.util import array_to_str, seqs_to_array,\
                                just_nucs

class TestJustNucs(TestCase):
    aln = LoadSeqs('data/just_nuc.fasta')
    
    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,
        seq1, seq6 to seq8 are expected to be eliminated by the code
        """
        d_aln = DenseAlignment(data=self.aln, MolType=DNA)
        
        nucs = just_nucs(d_aln.ArraySeqs)
        
        self.assertEqual(nucs,
                    [[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [2,0,1,2,3,1,2,2,2,0,2,2,2,2,2,3,0,2,2,3,0]])


class TestAlignSnpAnnotation(TestCase):
    d_aln = LoadSeqs('data/load_seqs_to_array.fasta', aligned=DenseAlignment, moltype=DNA)
    
    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,
        
        expect to convert seq0, seq2 to seq5 into numpy array
        """
        data = seqs_to_array(self.d_aln)
        expect = array([(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                        (1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                        (3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),
                        (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        (2,0,1,2,3,1,2,2,2,0,2,2,2,2,2,3,0,2,2,3,0)])
        
        self.assertEqual(data, expect)


class TestAlignSnpAnnotation(TestCase):
    seq_array = array([(2,0,1,2,2,1,2,0,2,0,2,2,2,2,2,3,3,2,2,2,0)])
    
    def test_array_to_str(self):
            """convert numpy array back to DNA sequence"""
            dna_str = array_to_str(self.seq_array)
            expect = ['ATCAACATATAAAAAGGAAAT']
            
            self.assertEqual(dna_str, expect)


if __name__ == '__main__':
    main()
