from numpy import fabs, array, around, floor, log10


import matplotlib.pyplot as plt
from pylab import gca
from cogent3 import DNA

import cogent3.core.moltype as moltype
from cogent3.draw.letters import letter_stack
from cogent3.util.union_dict import UnionDict

from mutation_motif.text import set_axis_scale, add_letter
from mutation_motif.util import FixedOrderFormatter

__author__ = "Jeremy Widman"
__copyright__ = "Copyright 2013, Jeremy Widman"
__credits__ = ["Jeremy Widman", "Yicheng Zhu", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


def est_ylim(char_heights):
    '''returns a ylim for character height axis plotting'''
    try:
        t = fabs(char_heights).sum(axis=1).max()
    except ValueError:
        t = fabs(char_heights).max()

    for i in range(1, 10):
        ylim = around(t, i)
        if ylim != 0:
            break

    if ylim < t:
        ylim *= 1.667
        ylim = around(ylim, i)

    ylim = max(ylim, 1e-6)

    return ylim


def draw_position(idx, idx_char_heights, sort_data=False, characters=None,
                  ax=None):
    """Draws characters at position.

        - idx is index in alignment
        - idx_char_heights is profile at index

    """
    if ax is None:
        ax = gca()

    if characters is None:
        characters = list(DNA)

    assert len(idx_char_heights) == len(characters)

    dx = 1.0
    x = idx
    y = 0
    data = list(zip(idx_char_heights, characters))
    if sort_data:
        data.sort()

    for i in range(len(characters)):
        dy, letter = data[i]
        add_letter(letter, x, y, dx, fabs(dy), ax=ax, invert=dy < 0)
        y += fabs(dy)

def set_anchored_ticks(ax, fontsize=14):
    """modifies the X-ticks to include negative numbers"""
    labels = ax.get_xticklabels()
    if not labels:
        return

    num = len(labels)
    mid = (num - 1) // 2
    assert 2 * mid + 1 == num, "Funny length"
    new_labels = [int(label.get_text()) - mid for label in labels]
    ax.set_xticklabels(new_labels, fontsize=fontsize)


def draw_alignment(char_heights, characters=None, ax=None, figsize=None,
                   ylim=None, figwidth=None, fig_callback=None,
                   set_ticks_func=set_anchored_ticks, fontsize=14,
                   verbose=False):
    """Takes in an alignment and creates an image from the alignment profile.
    """

    if ylim is None:
        # determine y-axis scale
        ylim = est_ylim(char_heights)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()

    fig = ax.get_figure()

    orig_len = len(char_heights)

    set_axis_scale(orig_len, ylim, ax=ax)

    if figwidth is None:
        figwidth = fig.get_figwidth()
        width_multiplier = orig_len / 20.
        fig.set_figwidth(width_multiplier * figwidth)
        print(fig.get_figwidth())
    else:
        fig.set_figwidth(figwidth)
        print(figwidth)

    ax.set_yticks(array([0.0, ylim / 2, ylim]))
    ax.set_yticklabels([0.0, ylim / 2, ylim], fontsize=14)
    ylabels = [ylabel.get_text() for ylabel in ax.get_yticklabels()]
    ax.set_yticklabels(ylabels, fontsize=14)
    ax.set_xticks(array([i + 0.5 for i in range(orig_len)]))
    ax.set_xticklabels([str(i) for i in range(orig_len)], rotation=-90)
    if set_ticks_func:
        set_ticks_func(ax, fontsize=fontsize)

    if fig_callback:
        fig_callback(fig)

    curr_idx = 0
    for i in range(orig_len):
        draw_position(i, char_heights[curr_idx],
                      characters=characters[i], ax=ax)
        curr_idx += 1

    return fig


def draw_multi_position(char_heights, characters, position_indices, ax=None,
                        figsize=None, ylim=None, figwidth=None,
                        fig_callback=None, set_ticks_func=set_anchored_ticks,
                        xtick_fontsize=14, ytick_fontsize=14, sort_data=False,
                        verbose=False):
    """Takes in an alignment and creates an image from the alignment profile.
    """
    if ylim is None:
        # determine y-axis scale
        ylim = est_ylim(char_heights)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()

    fig = ax.get_figure()

    orig_len = len(char_heights)

    set_axis_scale(orig_len, ylim, ax=ax)

    if figwidth is None:
        figwidth = fig.get_figwidth()
        width_multiplier = orig_len / 20.
        fig.set_figwidth(width_multiplier * figwidth)
    else:
        fig.set_figwidth(figwidth)

    ax.set_yticks(array([0.0, ylim / 2, ylim]))
    ax.set_yticklabels([0.0, ylim / 2, ylim])

    y_fmt = FixedOrderFormatter(floor(log10(ylim)))
    ax.yaxis.set_major_formatter(y_fmt)

    if fig_callback:
        fig_callback(fig)

    for position_index in position_indices:
        draw_position(position_index,
                      char_heights[position_index],
                      characters=characters[position_index], ax=ax,
                      sort_data=sort_data)

    ax.tick_params(axis='y', labelsize=ytick_fontsize, pad=ytick_fontsize // 2)
    ax.tick_params(axis='x', length=0)
    ax.set_xticks(array([i + 0.5 for i in range(orig_len)]))
    ax.set_xticklabels([str(i) for i in range(orig_len)], rotation=-90)
    if set_ticks_func:
        set_ticks_func(ax, fontsize=xtick_fontsize)

    return fig

def draw_multi_position_cogent3(char_heights, characters, position_indices, ax=None,
                                figsize=None, ylim=None, figwidth=None,
                                fig_callback=None, set_ticks_func=set_anchored_ticks,
                                xtick_fontsize=14, ytick_fontsize=14, sort_data=False,
                                verbose=False):
    """Takes in an alignment and creates an image from the alignment profile.
    """
    xaxis = "xaxis"+str(ax[1]+1 if ax[1] != 0 else "")
    yaxis = "yaxis"+str(ax[0]+1 if ax[0] != 0 else "")
    layout = UnionDict()

    if ylim is None:
        # determine y-axis scale
        ylim = est_ylim(char_heights)

    orig_len = len(char_heights)

    #set x and y range, set y ticks
    layout[xaxis] = dict(range=[0,orig_len], tickmode='array', tickvals=[i + 0.5 for i in range(orig_len)], ticktext=[str(i) for i in range(orig_len)], ticklen=0, tickangle=-90)
    layout[yaxis] = dict(range=[0,ylim], tickmode='array', tickvals=[0.0, ylim / 2, ylim], ticktext=[0.0, ylim / 2, ylim], tickfont=dict(size=ytick_fontsize))

    #y_fmt = FixedOrderFormatter(floor(log10(ylim)))
    #ax.yaxis.set_major_formatter(y_fmt)

    stack_data = [dict(zip(characters[i], char_heights[i])) for i in range(0, len(characters))]
    stacks = []
    for index, stack in enumerate(stack_data):
        middle, stack_shapes = letter_stack(stack, -2+index, 1, moltype.NT_COLORS, ax)
        stacks += stack_shapes

    layout["shapes"] = stacks

    #if set_ticks_func:
    #    set_ticks_func(ax, fontsize=xtick_fontsize)
    return layout
