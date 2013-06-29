from itertools import permutations
import os, glob

import numpy
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment
from mutation_motif import profile, util, height, logo, entropy, distribution

def PutN(data, direction):
    """put N text on plot"""
    N = data.shape[0]
    direction = r"$%s \rightarrow %s$" % tuple(direction.split('to'))
    def call(fig):
        ax = fig.gca()
        y = ax.get_ylim()[1]
        fig.gca().text(0.2, y * 0.75, direction, fontsize=20)
        fig.gca().text(0.2, y * 0.67, 'N=' + str(N))
    return call


fns = glob.glob('../data/snps_71/*.fasta.gz')
fns = [fn for fn in fns if 'to' in fn]

flank_size = 10
for fn in fns:
    direction = fn.split('-')[-1].split('.')[0]
    chosen_base = direction[0]
    basename = os.path.basename(fn).replace('.fasta.gz', '.png')
    outfile_name = os.path.join('../figs', basename)
    
    orig_seqs = util.load_from_fasta(fn)
    seqs = orig_seqs.ArraySeqs
    if seqs.shape[0] < 20:
        print 'Skipped, too few sequences', fn, seqs.shape
        continue
    
    orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base,
                                     step=3, flank_size=flank_size)

    ret = entropy.get_ret(orig, ctl)
    char_hts = height.get_re_char_heights(ret)
    callback = PutN(orig, direction)
    fig = logo.draw_alignment(char_hts.T, figsize=(9,3),
                            fig_callback=callback, verbose=False)
    if False:#orig.shape[0] < 5000:
        sigc = distribution.make_distribution(seqs, chosen_base, 3, flank_size, num_reps=100)
        fig.gca().plot(numpy.arange(0.5, sigc.shape[0], 1), sigc, color='k')
    
    fig.savefig(outfile_name, bbox_inches='tight')
    print 'Wrote', outfile_name
