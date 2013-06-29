from itertools import permutations
import os, glob

import numpy
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment
from mutation_motif import profile, util, height, logo, entropy, distribution

def PutN(data, direction, ylabel):
    """put N text on plot"""
    N = data.shape[0]
    direction = r"$%s \rightarrow %s$" % tuple(direction.split('to'))
    def call(fig):
        ax = fig.gca()
        ax.set_ylabel(ylabel, fontsize=20)
        y = ax.get_ylim()[1]
        fig.gca().text(0.2, y * 0.75, direction, fontsize=20)
        fig.gca().text(0.2, y * 0.67, 'N=' + str(N))
    return call

def get_re_char_hts(seqs, chosen_base, step, flank_size):
    orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base,
                                 step=step, flank_size=flank_size)

    ret = entropy.get_ret(orig, ctl)
    char_hts = height.get_re_char_heights(ret)
    return char_hts

def get_mi_char_hts(seqs, flank_size):
    """docstring for get_mi_char_hts"""
    midpt = (seqs.shape[1] -1) // 2
    assert 2 * midpt + 1 == seqs.shape[1], midpt
    start = midpt - flank_size
    seqs = seqs[:, start: start + (flank_size * 2 + 1)]
    probs = entropy.as_freq_matrix(seqs)
    mit = entropy.get_mit(seqs)
    mi = mit.sum(axis=0)
    char_hts = height.get_mi_char_heights(probs, mi)
    return char_hts

fns = glob.glob('../data/snps_71/rare-A-intro*.fasta.gz')
fns = [fn for fn in fns if 'to' in fn]
use_mi = True
ylabel = ['RE', 'MI'][use_mi]


flank_size = 5
for fn in fns:
    step = [1,3]['syn' in fn or 'missen' in fn]
    
    print 'working on', fn
    direction = fn.split('-')[-1].split('.')[0]
    chosen_base = direction[0]
    basename = os.path.basename(fn).replace('.fasta.gz', '.png')
    if use_mi:
        basename = 'mi-' + basename
    else:
        basename = 're-' + basename
        
    outfile_name = os.path.join('../figs', basename)
    
    orig_seqs = util.load_from_fasta(fn)
    seqs = orig_seqs.ArraySeqs
    seqs = util.just_nucs(seqs)
    if seqs.shape[0] < 200:
        print 'Skipped, too few sequences', fn, seqs.shape
        continue
    
    callback = PutN(seqs, direction, ylabel)
    if not use_mi:
        char_hts = get_re_char_hts(seqs, chosen_base, step, flank_size)
    else:
        char_hts = get_mi_char_hts(seqs, flank_size)
    
    fig = logo.draw_alignment(char_hts.T, figsize=(9,3),
                            fig_callback=callback, verbose=False)
    if seqs.shape[0] < 5000:
        sigc = distribution.make_distribution(seqs, chosen_base, 3, flank_size, num_reps=100)
        fig.gca().plot(numpy.arange(0.5, sigc.shape[0], 1), sigc, color='k')
    
    fig.savefig(outfile_name, bbox_inches='tight')
    print 'Wrote', outfile_name
