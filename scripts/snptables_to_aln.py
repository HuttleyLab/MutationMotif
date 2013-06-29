"""export seq files for different mutation types"""
from itertools import permutations

from mutation_motif.util import open_

directions = permutations('ACGT', 2)

is_autosome = lambda x: 'X' not in x and 'Y' not in x
is_xchrom = lambda x: 'X' in x

chrom = 'X'
freq = 0.1
freq_class = 'rare'

is_rare = lambda x: min(x) <= freq
is_common = lambda x: min(x) > freq

correct_freq = {'common': is_common, 'rare': is_rare}[freq_class]
correct_chrom = {'A': is_autosome, 'X': is_xchrom}[chrom]

for direction in directions:
    seen = set()
    effects = ['missense_variant', 'intron_variant', 'intergenic_variant', 'synonymous_variant']
    chroms = set()
    for effect in effects:
        name = (freq_class, chrom, effect,) + direction
        infilename = '../data/snps_71/%s_71.txt.gz' % effect
        outfilename = '../data/snps_71/%s-%s-%s_71-%sto%s.fasta.gz' % name
        
        with open_(infilename) as infile:
            with open_(outfilename, 'w') as outfile:
                num = 0
                for line in infile:
                    line = line.strip().split('\t')
                    coord = line[1]
                    
                    ancestor = line[6]
                    if len(ancestor) > 1:
                        continue
                    
                    alleles = line[5]
                    if ancestor not in alleles:
                        continue
                    
                    alleles = set(alleles.split('/'))
                    if len(alleles) != 2:
                        continue
                    
                    if not correct_chrom(coord):
                        continue
                    
                    got = alleles.difference(ancestor).pop()
                    if (ancestor, got) != direction:
                        continue
                    
                    label = line[0].strip()
                    if label in seen:
                        continue
                    
                    seen.update([label])
                    seq = line[7] + ancestor + line[8]
                    if 'N' in seq:
                        continue
                    
                    freqs = dict(eval(line[4]))
                    if not freqs:
                        continue
                    
                    if not correct_freq(freqs.values()):
                        continue
                    
                    record = '\n'.join(['>%s' % label, seq, ''])
                    outfile.write(record)
                    chroms.update([coord.split(':')[2]])
                    num += 1
                    if num % 1000 == 0:
                        print num
                    
                    if num == 1000000:
                        break
        
        print 'Wrote %s to %s' % (num, outfilename)
        print chroms
        print
