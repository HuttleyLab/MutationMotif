#!/usr/bin/env python
import os, time
from itertools import permutations, combinations
from optparse import make_option

import numpy
from matplotlib import pyplot

from cogent import LoadSeqs, DNA, LoadTable
from cogent.core.alignment import DenseAlignment
from cogent.util.option_parsing import parse_command_line_parameters
from cogent.maths.stats import chisqprob

from scitrack import CachingLogger

from mutation_motif import profile, util, logo, motif_count, log_lin

from mutation_motif.height import get_re_char_heights

LOGGER = CachingLogger(create_dir=True)

def make_summary(results):
    '''returns records from analyses as list'''
    rows = []
    for position_set in results:
        if type(position_set) != str:
            position = ':'.join(position_set)
        else:
            position = position_set
        
        re = results[position_set]['rel_entropy']
        dev = results[position_set]['deviance']
        df = results[position_set]['df']
        prob = results[position_set]['prob']
        formula = results[position_set]['formula']
        rows.append([position, re, dev, df, prob, formula])
    return rows

def get_selected_indices(stats):
    """returns indices for selecting dataframe records for display"""
    if 'group' in stats:
        indices = numpy.logical_and(stats['mut'] == 'M', stats['group'] == 1)
    else:
        indices = stats['mut'] == 'M'
    return indices

def get_grouped_combined_counts(table, position):
    """wraps motif_count.get_combined_counts for groups"""
    counts1 = motif_count.get_combined_counts(table.filtered("group in ('1', 1)"), position)
    counts2 = motif_count.get_combined_counts(table.filtered("group in ('2', 2)"), position)
    header = list(counts1.Header)
    counts1 = counts1.withNewColumn('group', lambda x : 1, columns=header[0])
    counts2 = counts2.withNewColumn('group', lambda x : 2, columns=header[0])
    header = ['group'] + header
    counts = LoadTable(header=header,
            rows=counts1.getRawData(header)+counts2.getRawData(header))
    counts.sorted(columns=['group', 'mut'])
    return counts

def get_position_effects(table, position_sets):
    pos_results = {}
    grouped = 'group' in table.Header
    if grouped:
        assert len(table.getDistinctValues('group')) == 2
    
    for position_set in position_sets:
        if not grouped:
            counts = motif_count.get_combined_counts(table, position_set)
        else:
            counts = get_grouped_combined_counts(table, position_set)
        
        rel_entropy, deviance, df, stats, formula = log_lin.position_effect(counts)
        if deviance < 0:
            p = 1.0
        else:
            p = chisqprob(deviance, df)
        
        pos_results[position_set] = dict(rel_entropy=rel_entropy,
                                        deviance=deviance, df=df, stats=stats,
                                        formula=formula, prob=p)
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
        mut_stats = stats[get_selected_indices(stats)][['base', 'ret']]
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
        mut_stats = stats[get_selected_indices(stats)][['base1', 'base2', 'ret']]
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
        mut_stats = stats[get_selected_indices(stats)][['base1', 'base2', 'base3', 'ret']]
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
    mut_stats = stats[get_selected_indices(stats)][['base1', 'base2', 'base3', 'base4', 'ret']]
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
     make_option('-1','--countsfile', help='tab delimited file of counts.'),
     make_option('-o','--outpath', help='Directory path to write data.'),
    ]

script_info['optional_options'] = [
    make_option('-2','--countsfile2',
        help='second group motif counts file.'),
    make_option('--format', default='pdf', choices=['pdf', 'png'],
        help='Plot format.'),
    make_option('-F','--force_overwrite', action='store_true', default=False,
        help='Overwrite existing files.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'

def single_group(opts):
    option_parser, opts, args =\
       parse_command_line_parameters(disallow_positional_arguments=False, **script_info)
    
    outpath = util.abspath(opts.outpath)
    if not opts.dry_run:
        util.create_path(outpath)
        runlog_path = os.path.join(outpath, "analysis.log")
        LOGGER.log_file_path = runlog_path
        LOGGER.write(str(vars(opts)), label='vars')
    
    counts_filename = util.abspath(opts.countsfile)
    counts_table = LoadTable(counts_filename, sep='\t')
    
    
    LOGGER.input_file(counts_filename, label="countsfile1_path")
    
    num_pos = sum(1 for c in counts_table.Header if c.startswith('pos'))
    if num_pos != 4:
        raise ValueError("Requires four positions for analysis")
    
    if opts.countsfile2:
        print "Performing 2 group analysis"
        counts_table1 = counts_table.withNewColumn('group', lambda x: '1',
                                    columns=counts_table.Header[0])
        
        fn2 = util.abspath(opts.countsfile2)
        counts_table2 = LoadTable(fn2, sep='\t')
        
        LOGGER.input_file(fn2, label="countsfile2_path")
        
        counts_table2 = counts_table2.withNewColumn('group', lambda x: '2',
                                    columns=counts_table2.Header[0])
        # now combine
        header = ['group'] + counts_table2.Header[:-1]
        raw1 = counts_table1.getRawData(header)
        raw2 = counts_table2.getRawData(header)
        counts_table = LoadTable(header=header, rows=raw1+raw2)
        
        if not opts.dry_run:
            outfile = os.path.join(outpath, 'group_counts_table.txt')
            counts_table.writeToFile(outfile, sep='\t')
            LOGGER.output_file(outfile, label="group_counts")
            
    
    if opts.dry_run or opts.verbose:
        print
        print counts_table
        print
    
    positions = [c for c in counts_table.Header if c.startswith('pos')]
    
    # Collect statistical analysis results
    summary = []
    
    max_results = {}
    # Single position analysis
    print "Doing single position analysis"
    single_results = single_position_effects(counts_table, positions)
    summary += make_summary(single_results)
    
    max_results[1] = max(single_results[p]['rel_entropy'] for p in single_results)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "1.json")
        util.dump_loglin_stats(single_results, outfilename)
        LOGGER.output_file(outfilename, label="analysis1")
    
    fig = get_single_position_fig(single_results, positions, (9,3), fig_width=2.25)
    fig.tight_layout()
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "1.pdf")
        fig.savefig(outfilename)
        print "Wrote", outfilename
        fig.clf() # refresh for next section
    
    print "Doing two positions analysis"
    results = get_two_position_effects(counts_table, positions)
    summary += make_summary(results)
    
    max_results[2] = max(results[p]['rel_entropy'] for p in results)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "2.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis2")
    
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
    summary += make_summary(results)
    
    max_results[3] = max(results[p]['rel_entropy'] for p in results)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "3.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis3")

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
    summary += make_summary(results)
    
    max_results[4] = max(results[p]['rel_entropy'] for p in results)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "4.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis4")
    
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
    
    summary = LoadTable(header=['Position', 'RE', 'Deviance', 'df', 'prob', 'formula'],
            rows=summary)
    if not opts.dry_run:
        outfilename = os.path.join(outpath, "summary.txt")
        summary.writeToFile(outfilename, sep='\t')
        LOGGER.output_file(outfilename, label="summary")
    
    print summary
    
    print "Done! Check %s for your results" % outpath
    

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(disallow_positional_arguments=False, **script_info)
    
    single_group(opts)
