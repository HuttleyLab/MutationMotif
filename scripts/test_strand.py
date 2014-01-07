import sys, os

from cogent.util.unit_test import TestCase, main
from strand import reverse_complement_record, get_rc_record
from cogent.db.ensembl import HostAccount, Genome

account = HostAccount('gduserv.anu.edu.au', 'compgen', 'compgen')
human = Genome('human', Release=71, account=account)


class TestStrand(TestCase):
    def est_reverse_record(self):
        """correctly identify when a record needs to be reversed"""
        # gene is on plus strand and snp is on plus strand
        self.assertFalse(reverse_complement_record(gene_strand=1, snp_strand=1))
        # gene is on plus strand and snp is on minus strand
        self.assertTrue(reverse_complement_record(gene_strand=1, snp_strand=-1))
        # gene is on minus strand and snp is on plus strand
        self.assertFalse(reverse_complement_record(gene_strand=-1, snp_strand=1))
        # gene is on minus strand and snp is on minus strand
        self.assertTrue(reverse_complement_record(gene_strand=-1, snp_strand=-1))
        # invalid data
        self.assertRaises(RuntimeError, reverse_complement_record, 0, 0)
    
    def test_flanks_are_neg_strand(self):
        """establish relationship between strand of SNP and gene"""
        # SNPs reflecting different possible relationships
        # strand of gene_snp
        plus_plus='rs2462492'
        plus_minus='rs1063798'
        minus_minus='rs3103691'
        minus_plus='rs200394590'
        
        # flank is in same orientation as gene
        snp = list(human.getVariation(Symbol=plus_plus))[0] # don't rc for flank
        flank5, flank3 = snp.FlankingSeq
        gene = list(snp.getFeatures('gene'))[0]
        self.assertTrue(str(flank5) in str(gene.Seq))
        
        snp = list(human.getVariation(Symbol=minus_minus))[0] # rc for flank
        flank5, flank3 = snp.FlankingSeq
        gene = list(snp.getFeatures('gene'))[0]
        self.assertTrue(str(flank5) in str(gene.Seq))
        
        # flank is in different orientation to gene
        snp = list(human.getVariation(Symbol=plus_minus))[0] # rc for flank
        flank5, flank3 = snp.FlankingSeq
        gene = list(snp.getFeatures('gene'))[0]
        self.assertTrue(str(flank5) in str(gene.Seq.rc()))
        
        snp = list(human.getVariation(Symbol=minus_plus))[0] # rc for flank
        flank5, flank3 = snp.FlankingSeq
        gene = list(snp.getFeatures('gene'))[0]
        self.assertTrue(str(flank5) in str(gene.Seq.rc()))
    
    def test_rc_record(self):
        """reverse complement a record"""
        alleles = set('TC')
        ancestor = 'T'
        af = {'C': 0.9885, 'T': 0.0115}
        f_seq_5 = 'ATAAGATCTTATGTTCTT'
        f_seq_3 = 'GACATTGCAGGAAAATCT'
        alleles_rc, ancestor_rc, allele_freqs_rc, flank_5_rc, flank_3_rc = \
                    get_rc_record(alleles, ancestor, af, f_seq_5, f_seq_3)
        self.assertEqual(alleles_rc, set('AG'))
        self.assertEqual(ancestor_rc, 'A')
        expect_freqs = {'A': 0.0115, 'G': 0.9885}
        self.assertEqual(allele_freqs_rc, expect_freqs)
    

if __name__ == "__main__":
    main()
