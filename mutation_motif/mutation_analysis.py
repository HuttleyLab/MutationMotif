#!/usr/bin/env python
import os, time
from itertools import permutations, combinations
from optparse import make_option

import numpy
from matplotlib import pyplot

from cogent import LoadSeqs, DNA, LoadTable
from cogent.core.alignment import DenseAlignment
from cogent.util.option_parsing import parse_command_line_parameters

from mutation_motif import profile, util, logo, motif_count, log_lin

from mutation_motif.height import get_re_char_heights

def get_position_effects(table, position_sets):
    pos_results = {}
    for position_set in position_sets:
        counts = motif_count.get_combined_counts(table, position_set)
        rel_entropy, deviance, df, stats = log_lin.position_effect(counts)
        pos_results[position_set] = dict(rel_entropy=rel_entropy,
                                        deviance=deviance, df=df, stats=stats)
    return pos_results

def single_position_effects(table, positions):
    single_results = get_position_effects(table, positions)
    return single_results

def get_single_position_fig(single_results, positions, figsize, fig_width=None, fontsize=14):
    num_pos = len(positions) + 1
    mid = num_pos // 2

    position_re = numpy.zeros((num_pos,), float)
    rets = numpy.zeros((4, num_pos), float)
    characters = [list('ACGT') for i in range(num_pos)]
    for index, pos in enumerate(positions):
        if index >= mid:
            index += 1

        stats = single_results[pos]['stats']
        position_re[index] = single_results[pos]['rel_entropy']
        mut_stats = stats[stats['mut'] == 'M'][['base', 'ret']]
        mut_stats = mut_stats.sort('ret')
        characters[index] = list(mut_stats['base'])
        rets[:,index] = mut_stats['ret']

    heights = get_re_char_heights(rets, re_positionwise=position_re)
    fig = logo.draw_multi_position(heights.T, characters=characters, position_indices=range(num_pos),
                               figsize=figsize, fig_width=fig_width, verbose=False)
    
    if fig_width:
        fig.set_figwidth(fig_width)
    
    ax = fig.gca()
    ax.set_xlabel('Position', fontsize=fontsize)
    ax.set_ylabel('RE', rotation='vertical', fontsize=fontsize)
    return fig

def get_resized_array_coordinates2(positions, motifs):
    '''coords for the position combinations for the smaller array used for plotting'''
    # to avoid blank rows/columns in the trellis plot, we reduce the array size
    num_pos = len(positions)
    a = numpy.zeros((num_pos, num_pos), object)
    for motif in motifs:
        indices = map(positions.index, motif)
        indices.reverse()
        a[indices[0]][indices[1]] = motif
    
    new = a[(a != 0).any(axis=1)]
    # it's a 2D array now
    new = new[:, (new != 0).any(axis=0)]
    mapping = {}
    for i in range(new.shape[0]):
        for j in range(new.shape[1]):
            if new[i, j] != 0:
                mapping[new[i, j]] = (i, j)
    return mapping

def get_two_position_effects(table, positions):
    two_pos_results = get_position_effects(table, list(combinations(positions, 2)))
    return two_pos_results

def get_two_position_fig(two_pos_results, positions, figsize, fig_width=None, fontsize=14):
    position_sets = list(combinations(positions, 2))
    array_coords = get_resized_array_coordinates2(positions, position_sets)
    coords = array_coords.values()
    xdim = max(v[0] for v in coords) + 1
    ydim = max(v[1] for v in coords) + 1
    fig, axarr = pyplot.subplots(xdim, ydim, figsize=figsize, sharex=True, sharey=True)

    for i in range(xdim-1):
        for j in range(i+1, ydim):
            ax = axarr[i,j]
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    rel_entropies = []
    for position_set in position_sets:
        rel_entropies.append(two_pos_results[position_set]['rel_entropy'])
    ylim = logo.est_ylim(rel_entropies)
    
    num_pos = len(positions) + 1

    mid_pos = num_pos // 2

    position_re = numpy.zeros((num_pos,), float)
    multi_positions = {}
    characters = numpy.zeros((num_pos,  16), str)

    for pair in position_sets:
        rets = numpy.zeros((16, num_pos), float)
        indices = map(positions.index, pair)
        row, col = array_coords[pair]
        ax = axarr[row, col]

        # now adjust indices to reflect position along sequence
        for i in range(2):
            if indices[i] >= mid_pos:
                indices[i] += 1

        position_re.put(indices, two_pos_results[pair]['rel_entropy'])

        stats = two_pos_results[pair]['stats']
        mut_stats = stats[stats['mut'] == 'M'][['base1', 'base2', 'ret']]
        mut_stats = mut_stats.sort(columns='ret')

        characters[indices[0]] = list(mut_stats['base1'])
        characters[indices[1]] = list(mut_stats['base2'])

        for index in indices:
            rets[:, index] = mut_stats['ret']

        heights = get_re_char_heights(rets, re_positionwise=position_re)
        multi_positions[pair] = dict(rets=rets, indices=indices, characters=characters, heights=heights)
        logo.draw_multi_position(char_heights=heights.T, characters=characters,
                                 position_indices=indices, ax=ax, ylim=ylim)
    return fig

def get_resized_array_coordinates3(positions, position_set):
    '''coords for the position combinations for the smaller array used for plotting'''
    # to avoid blank rows/columns in the trellis plot, we reduce the array size
    num_pos = len(positions)
    a = numpy.zeros((num_pos, num_pos, num_pos), object)
    for triple in position_set:
        indices = map(positions.index, triple)
        indices.reverse()
        a[indices[0]][indices[1]][indices[2]] = triple
    
    new = a.flatten()
    new = new[new != 0]
    new.sort()
    new = new.reshape((2,2))
    mapping = {}
    for i in range(new.shape[0]):
        for j in range(new.shape[1]):
            if new[i, j] != 0:
                mapping[new[i, j]] = (i, j)
    return mapping

def get_three_position_effects(table, positions):
    three_pos_results = get_position_effects(table, list(combinations(positions, 3)))
    return three_pos_results

def get_three_position_fig(three_pos_results, positions, figsize, fig_width=None, fontsize=14):
    position_sets = list(combinations(positions, 3))
    array_coords = get_resized_array_coordinates3(positions, position_sets)

    coords = array_coords.values()
    xdim = max(v[0] for v in coords) + 1
    ydim = max(v[1] for v in coords) + 1

    fig, axarr = pyplot.subplots(xdim, ydim, figsize=(9,9), sharex=True, sharey=True)

    for i in range(xdim):
        for j in range(ydim):
            if (i,j) in coords:
                continue

            ax = axarr[i,j]
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)


    num_pos = len(positions) + 1
    mid_pos = num_pos // 2
    
    rel_entropies = []
    for position_set in position_sets:
        rel_entropies.append(three_pos_results[position_set]['rel_entropy'])
    ylim = logo.est_ylim(rel_entropies)



    position_re = numpy.zeros((num_pos,), float)
    multi_positions = {}
    characters = numpy.zeros((num_pos, 64), str)

    for motif in combinations(positions, 3):
        rets = numpy.zeros((64, num_pos), float)
        indices = map(positions.index, motif)
        row, col = array_coords[motif]
        ax = axarr[row, col]

        # now adjust indices to reflect position along sequence
        for i in range(len(indices)):
            if indices[i] >= mid_pos:
                indices[i] += 1

        position_re.put(indices, three_pos_results[motif]['rel_entropy'])

        stats = three_pos_results[motif]['stats']
        mut_stats = stats[stats['mut'] == 'M'][['base1', 'base2', 'base3', 'ret']]
        mut_stats = mut_stats.sort(columns='ret')

        characters[indices[0]] = list(mut_stats['base1'])
        characters[indices[1]] = list(mut_stats['base2'])
        characters[indices[2]] = list(mut_stats['base3'])


        for index in indices:
            rets[:, index] = mut_stats['ret']

        heights = get_re_char_heights(rets, re_positionwise=position_re)
        multi_positions[motif] = dict(rets=rets, indices=indices, characters=characters, heights=heights)
        logo.draw_multi_position(char_heights=heights.T, characters=characters,
                                 position_indices=indices, ax=ax, ylim=ylim)

    return fig

def get_four_position_effects(table, positions):
    result = get_position_effects(table, list(combinations(positions, 4)))
    return result

def get_four_position_fig(four_pos_results, positions, figsize, fig_width=None, fontsize=14):
    position_sets = list(combinations(positions, 4))
    assert len(position_sets) == 1
    rel_entropies = [four_pos_results[position_sets[0]]['rel_entropy']]
    ylim = logo.est_ylim(rel_entropies)
    
    rel_entropy = rel_entropies[0]
    
    fig = pyplot.figure(figsize=(9,9))
    ax = fig.gca()

    num_pos = len(positions) + 1

    mid_pos = num_pos // 2

    position_re = numpy.zeros((num_pos,), float)
    characters = numpy.zeros((num_pos, 256), str)

    rets = numpy.zeros((256, num_pos), float)
    indices = range(4)

    # now adjust indices to reflect position along sequence
    for i in range(len(indices)):
        if indices[i] >= mid_pos:
            indices[i] += 1

    position_re.put(indices, rel_entropy)
    stats = four_pos_results[position_sets[0]]['stats']
    mut_stats = stats[stats['mut'] == 'M'][['base1', 'base2', 'base3', 'base4', 'ret']]
    mut_stats = mut_stats.sort(columns='ret')

    characters[indices[0]] = list(mut_stats['base1'])
    characters[indices[1]] = list(mut_stats['base2'])
    characters[indices[2]] = list(mut_stats['base3'])
    characters[indices[3]] = list(mut_stats['base4'])

    for index in indices:
        rets[:, index] = mut_stats['ret']

    heights = get_re_char_heights(rets, re_positionwise=position_re)
    logo.draw_multi_position(char_heights=heights.T, characters=characters,
                             position_indices=indices, ax=ax, ylim=ylim)

    return fig

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = \
"log-linear analysis of neighbouring bases."\
+ " Writes counts tables, estimated statistics and figures to the specified"\
+ "directory outpath."

script_info['required_options'] = [
     make_option('-a','--alignfile', help='fasta aligned file centred on mutated position.'),
     make_option('-o','--outpath', help='Directory path to write data.'),
    make_option('--direction', default=None,
     choices=['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT', 'GtoA', 'GtoC',
             'GtoT', 'TtoA', 'TtoC', 'TtoG'], help='Mutation direction.'),
    ]

script_info['optional_options'] = [
    make_option('-S', '--seed',
        help='Seed for random number generator (e.g. 17, or 2015-02-13). Defaults to system time.'),
    make_option('--coding', action='store_true', default=False,
        help='If True, samples pseudo-SNPs at 3-bp intervals from real SNP.'),
    make_option('--format', default='pdf', choices=['pdf', 'png'],
        help='Plot format.'),
    make_option('-F','--force_overwrite', action='store_true', default=False,
        help='Overwrite existing files.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(disallow_positional_arguments=False, **script_info)
    
    flank_size = 2
    
    outpath = util.abspath(opts.outpath)
    if not opts.dry_run:
        util.create_path(outpath)
    
    counts_filename = os.path.join(outpath, "counts_table.txt")
    
    if not opts.seed:
        opts.seed = str(time.time())
        print "NOTE: set random number seed to '%s'" % (opts.seed)
    
    if not os.path.exists(counts_filename):
        print "Deriving counts from sequence file"
        direction = tuple(opts.direction.split('to'))
        chosen_base = direction[0]
    
        fn = util.abspath(opts.alignfile)
        orig_seqs = util.load_from_fasta(fn)
        seqs = orig_seqs.ArraySeqs
        seqs = util.just_nucs(seqs)
        orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base, step=1,
                                        flank_size=flank_size, seed=opts.seed)
    
        # convert profiles to a motif count table
        orig_counts = motif_count.profile_to_seq_counts(orig, flank_size=flank_size)
        ctl_counts = motif_count.profile_to_seq_counts(ctl, flank_size=flank_size)
        counts_table = motif_count.get_count_table(orig_counts, ctl_counts, flank_size*2)
        counts_table = counts_table.sorted(columns='mut')
        if not opts.dry_run:
            counts_table.writeToFile(counts_filename, sep='\t')
    else:
        print "Loading existing counts_table"
        counts_table = LoadTable(counts_filename, sep='\t')
    
    
    if opts.dry_run or opts.verbose:
        print
        print counts_table
        print
    
    positions = [c for c in counts_table.Header if c.startswith('pos')]
    
    max_results = {}
    # Single position analysis
    print "Doing single position analysis"
    single_results = single_position_effects(counts_table, positions)
    max_results[1] = max(single_results[p]['rel_entropy'] for p in single_results)
    fig = get_single_position_fig(single_results, positions, (9,3), fig_width=2.25)
    fig.tight_layout()
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "1.pdf")
        fig.savefig(outfilename)
        print "Wrote", outfilename
        fig.clf() # refresh for next section
    
    print "Doing two positions analysis"
    results = get_two_position_effects(counts_table, positions)
    max_results[2] = max(results[p]['rel_entropy'] for p in results)
    fig = get_two_position_fig(results, positions, (9,9))
    fig.set_figwidth(9)
    fig.text(0.5, 0.05, 'Position', ha='center', va='center', fontsize=14)
    fig.text(0.01, 0.5, 'RE', ha='center', va='center', rotation='vertical', fontsize=14)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "2.pdf")
        fig.savefig(outfilename)
        print "Wrote", outfilename
        fig.clf() # refresh for next section
    
    print "Doing three positions analysis"
    results = get_three_position_effects(counts_table, positions)
    max_results[3] = max(results[p]['rel_entropy'] for p in results)
    fig = get_three_position_fig(results, positions, (9,9))
    fig.set_figwidth(9)
    fig.text(0.5, 0.05, 'Position', ha='center', va='center', fontsize=14)
    fig.text(0.01, 0.5, 'RE', ha='center', va='center', rotation='vertical', fontsize=14)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "3.pdf")
        fig.savefig(outfilename)
        print "Wrote", outfilename
        fig.clf() # refresh for next section
    
    print "Doing four positions analysis"
    results = get_four_position_effects(counts_table, positions)
    max_results[4] = max(results[p]['rel_entropy'] for p in results)
    fig = get_four_position_fig(results, positions, (9,9))
    fig.set_figwidth(9)
    ax = fig.gca()
    ax.set_xlabel('Position', fontsize=14)
    ax.set_ylabel('RE', fontsize=14)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "4.pdf")
        fig.savefig(outfilename)
        print "Wrote", outfilename
        fig.clf() # refresh for next section
    
    # now generate summary plot
    bar_width = 0.5
    index = numpy.arange(4)
    bar = pyplot.bar(index, [max_results[i] for i in range(1,5)], bar_width)
    pyplot.xticks(index+(bar_width/2.), range(1,5))
    pyplot.ylabel("RE$_{max}$")
    pyplot.xlabel("Motif Length")
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "summary.pdf")
        pyplot.savefig(outfilename)
        print "Wrote", outfilename
    
    print "Done! Check %s for your results" % outpath