from numpy import matrix, fabs, array, ceil, around

from matplotlib import pyplot
pyplot.switch_backend('Agg')

import matplotlib.pyplot as plt
from pylab import savefig, gcf, gca, clf
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText

from cogent.align.weights.util import AlnToProfile, DNA, DNA_ORDER

from mutation_motif.text import set_axis_scale, add_letter

__author__ = "Jeremy Widman"
__copyright__ = "Copyright 2013, Jeremy Widman"
__credits__ = ["Jeremy Widman", "Yicheng Zhu", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def draw_position(idx, idx_char_heights):
    """Draws characters at position.
        
        - idx is index in alignment
        - idx_char_heights is profile at index
        
    """
    dx = 1.0
    x = idx
    y = 0
    data = zip(idx_char_heights, DNA)
    data.sort()
    for i in range(4):
        dy, letter = data[i]
        add_letter(letter, x, y, dx, fabs(dy), invert=dy<0)
        y += fabs(dy)

def set_anchored_ticks(fig):
    """modifies the X-ticks to include negative numbers"""
    ax = fig.gca()
    labels = ax.get_xticklabels()
    num = len(labels)
    mid = (num - 1) / 2
    assert 2 * mid + 1 == num, "Funny length"
    new_labels = [int(label.get_text()) - mid for label in labels]
    ax.set_xticklabels(new_labels, fontsize=20)

def draw_alignment(char_heights, ax=None, figsize=None, ylim=None, fig_callback=None, set_ticks_func=set_anchored_ticks, verbose=False):
    """Takes in an alignment and creates an image from the alignment profile.
    """
    
    if ylim is None:
        # determine y-axis scale
        t = fabs(char_heights).sum(axis=1).max()
        for i in range(1, 10):
            ylim = around(t, i)
            if ylim != 0:
                break
        
        if verbose:
            print t, ylim
        
        if ylim < t:
            ylim *= 1.667
            ylim = around(ylim, i)
    
    mid = (char_heights.shape[1] -1 ) / 2
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
    
    fig = ax.get_figure()
    
    orig_len = len(char_heights)
    
    set_axis_scale(orig_len, ylim)
    
    fig_width = fig.get_figwidth()
    width_multiplier = orig_len / 20.
    fig.set_figwidth(width_multiplier*fig_width)
    
    ax.set_yticks(array([0.0, ylim / 2, ylim]))
    ax.set_yticklabels([0.0, ylim / 2, ylim], fontsize=20)
    ylabels = [ylabel.get_text() for ylabel in ax.get_yticklabels()]
    ax.set_yticklabels(ylabels, fontsize=20)
    ax.set_xticks(array([i+0.5 for i in range(orig_len)]))
    ax.set_xticklabels([str(i) for i in range(orig_len)], rotation=-90)
    if set_ticks_func:
        set_ticks_func(fig)
    
    if fig_callback:
        fig_callback(fig)
    
    curr_idx = 0
    for i in range(orig_len):
        draw_position(i, char_heights[curr_idx])
        curr_idx += 1
    
    gca().yaxis.grid(True)
    return fig

