from cogent import LoadTable, DNA

def reverse_complement_record(gene_strand, snp_strand):
    """returns True if the SNP record needs to be reverse complemented to put on transcribed strand"""
    
    if gene_strand == snp_strand:
        # gene is + snp is +; False
        rc = False
    elif gene_strand != snp_strand:
        # gene is - snp is -; True
        rc = False
    
    return rc

def get_rc_record(alleles, ancestor, allele_freqs, flank_5, flank_3):
    """reverse complements the alleles, ancestror, flanking seqs, and allele freqs
    """
    complement = DNA.complement
    alleles_rc = set([complement(b) for b in alleles])
    ancestor_rc = complement(ancestor)
    
    allele_freqs_rc = {}
    for allele, freq in allele_freqs.items():
        allele_freqs_rc[complement(allele)] = freq
    
    flank_5_rc = str(DNA.makeSequence(flank_5).rc())
    flank_3_rc = str(DNA.makeSequence(flank_3).rc())
    return alleles_rc, ancestor_rc, allele_freqs_rc, flank_5_rc, flank_3_rc
