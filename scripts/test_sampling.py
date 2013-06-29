import os, sys

sys.path.append('..')

from cogent import DNA
from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl import Species, HostAccount, Genome

from mutation_motif import sample_snp

account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
human = Genome('human', Release=71, account=account)

class TestSnpSampling(TestCase):
    def est_validated_only(self):
        """correctly discriminates validated, non-validated SNPs"""
        symbol = 'rs2241629' # should fail
        snps = list(human.getVariation(Symbol=symbol))[0]
        print snps.AlleleFreqs
        sample_snp.get_max_allele_freqs(snps.AlleleFreqs)
    
    def test_snp_variant_classification(self):
        """correctly classify SNPs"""
        pass
    
    def test_snp_variant_freq(self):
        """correctly extract minimum allele frequency"""
        pass
    
    def test_ancestral(self):
        """correctly idenfity mutation direction"""
        pass
    
    def test_snp_strand(self):
        """correctly establish SNP strand for the different classes of SNPs"""
        pass
    
    def test_coorect_coords(self):
        """SNPs have correct genomic coordinates"""
        pass
    
    def test_infer_mutation_direction(self):
        """correctly identify the mutation direction"""
        pass
    
    def est_snp_data_dump(self):
        """returns correct data for saving, handles reverse complement"""
        # [snp.Symbol, rc, str(snp.Location), str(snp.Location.Strand),
        #             snp.Effect, str(allele_freqs), alleles, str(ancestral),
        #             str(flank_5prime), str(flank_3prime)]
        # SNPs are for BRCA1, ENSG00000012048, which is on -ve strand
        # rs6554830 is a syn A/G SNP
        # 5' Near Seq 30 bp: 
        # 3' Near Seq 30 bp: 
        snp = list(human.getVariation(Symbol='rs6554830'))[0]
        print '+ strand'
        got = sample_snp.get_snp_dump_data(snp) # forward strand
        print got
        print snp
        print snp.Ancestral
        # alleles
        self.assertEqual(set(got[6].split('/')), set(['A', 'G']))
        # ancestral allele
        self.assertEqual(got[7], 'G')
        flanks = DNA.makeSequence(got[8]+got[9])
        print '- strand'
        got = sample_snp.get_snp_dump_data(snp, rc=True) # forward strand
        print got
        self.assertEqual(set(got[6].split('/')), set(['T', 'C']))
        self.assertEqual(got[8]+got[9], str(flanks.rc()))
        self.assertEqual(got[7], 'C')
        
        
        # rs3810870 is a syn G/C SNP
        # 5' Near Seq 30 bp: CGCTCGTCCAGGGAGGCCGGACAGTGGTTG
        # 3' Near Seq 30 bp: GCAGGTTCCCAGGCACCACCCCCACTGTCC
        snp = list(human.getVariation(Symbol='rs3810870'))[0]
        print '+ strand'
        got = sample_snp.get_snp_dump_data(snp) # forward strand
        print got
        # alleles
        self.assertEqual(set(got[6].split('/')), set(['C', 'G']))
        # ancestral allele -- based on inspection, but seems to differ to dbSNP
        self.assertEqual(got[7], 'C')
        flanks = DNA.makeSequence(got[8]+got[9])
        print '- strand'
        got = sample_snp.get_snp_dump_data(snp, rc=True) # forward strand
        print got
        self.assertEqual(set(got[6].split('/')), set(['C', 'G']))
        self.assertEqual(got[8]+got[9], str(flanks.rc()))
        self.assertEqual(got[7], 'G')
        

if __name__ == '__main__':
    main()
