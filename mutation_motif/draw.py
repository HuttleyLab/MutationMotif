import os
import sys
from configparser import SafeConfigParser

import click
import numpy
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MultipleLocator
from cogent3 import DNA
from scitrack import CachingLogger

from mutation_motif import util, logo, text, entropy
from mutation_motif.height import get_re_char_heights, get_mi_char_heights

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


LOGGER = CachingLogger(create_dir=True)


def format_float(cutoff, float_places=3, sci_places=1):
    """factory function for swapping from float to scientific notation"""
    float_template = "%%.%df" % float_places
    sci_template = "%%.%de" % sci_places

    def call(value, *args, **kwargs):
        if abs(value) == 0:
            txt = float_template % value
        elif abs(value) < cutoff:
            txt = sci_template % value
        else:
            txt = float_template % value
        return txt

    return call


def get_selected_indices(stats):  # from mutation_analysis
    """returns indices for selecting dataframe records for display"""
    if 'group' in stats:
        indices = numpy.logical_and(stats['mut'] == 'M', stats['group'] == 1)
    else:
        indices = stats['mut'] == 'M'
    return indices


def get_plot_data(single_results, positions):
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
        mut_stats = mut_stats.sort_values(by='ret')
        characters[index] = list(mut_stats['base'])
        rets[:, index] = mut_stats['ret']

    heights = get_re_char_heights(rets, re_positionwise=position_re)
    return heights.T, characters, list(range(num_pos))


def draw_position_grid(directions, sample_size=False, width=8, height=8,
                       title_space=1.1, axis_font_size=20, tick_font_size=10,
                       ylim=None):
    """docstring for draw_position_grid"""
    f, axes = pyplot.subplots(4, 4, sharex=True, sharey=True,
                              figsize=(width, height))

    bases = list('CTAG')
    positions = None
    plottables = []
    adaptive_y = 0
    for direction in directions:
        data = directions[direction]
        if positions is None:
            positions = list(data.keys())
            positions.sort()
        number = data[positions[0]]['stats']["count"].sum() / 2
        heights, characters, indices = get_plot_data(data, positions)
        adaptive_y = max(adaptive_y, logo.est_ylim(heights))
        plottables.append([direction, heights, characters, indices, number])

    if ylim is None:
        ylim = adaptive_y

    for direction, heights, characters, indices, number in plottables:
        fr, to = list(map(bases.index, direction.split('to')))
        ax = axes[fr, to]
        fig = logo.draw_multi_position(heights, characters=characters,
                                       position_indices=indices, ylim=ylim,
                                       ax=ax, figwidth=width, verbose=False)
        if sample_size:
            y = ax.get_ylim()[1]
            ax.text(0.2, y * 0.85, "N={:,}".format(number), fontsize=10)

    xformat = FuncFormatter(format_float(1e-3, float_places=2))

    for i in range(4):
        top_ax = axes[0, i]
        top_ax.set_title(
            bases[i], fontsize=axis_font_size, weight="bold", y=1.1)

        lft_ax = axes[i, 0]
        for yticklabel in lft_ax.get_yticklabels():
            yticklabel.set_fontsize(tick_font_size)
            yticklabel.set_rotation(0)
        lft_ax.yaxis.set_major_formatter(FuncFormatter(xformat))
        lft_ax.set_ylabel(bases[i], rotation=0,
                          fontsize=axis_font_size, weight="bold")
        lft_ax.yaxis.labelpad = axis_font_size

        btm_ax = axes[-1, i]
        for xticklabel in btm_ax.get_xticklabels():
            xticklabel.set_fontsize(tick_font_size)
            xticklabel.set_rotation(0)

    f.tight_layout()

    return f


_figpath = click.option(
    '--figpath', help='Filename for plot file. Overides format.')
_plot_cfg = click.option('--plot_cfg',
                         help='Config file for plot size, font size settings.')
_format = click.option('--format', default='pdf',
                       type=click.Choice(['pdf', 'png']),
                       help='Plot figure format.')
_sample_size = click.option('--sample_size', is_flag=True,
                            help='Include sample size on each subplot.')
_no_type3 = util.no_type3_font
_force_overwrite = click.option('-F', '--force_overwrite', is_flag=True,
                                help='Overwrite existing files.')
_dry_run = click.option('-D', '--dry_run', is_flag=True,
                        help='Do a dry run of the analysis without writing '
                        'output.')


@click.group()
def main():
    """draw mutation motif logo's and spectra"""
    pass


_paths_cfg = click.option('--paths_cfg', required=True,
                          help='Text file listing path for 1.json file for '
                          'each mutation direction (e.g. AtoG).')


@main.command()
@_paths_cfg
@_plot_cfg
@_figpath
@_format
@_no_type3
@_sample_size
@_force_overwrite
@_dry_run
def nbr_grid(paths_cfg, plot_cfg, figpath, format, no_type3, sample_size,
             force_overwrite, dry_run):
    '''draws grid of sequence logo's from neighbour analysis'''
    if no_type3:
        util.exclude_type3_fonts()

    args = locals()
    LOGGER.log_message(str(args), label='vars')

    config_path = util.abspath(paths_cfg)
    indir = os.path.dirname(config_path)
    parser = SafeConfigParser()
    parser.optionxform = str  # stops automatic conversion to lower case
    parser.read(config_path)

    json_paths = {}
    for direction, path in parser.items("json_paths"):
        # assumes paths are relative to indir
        path = os.path.join(indir, path)
        if not os.path.exists(path):
            print("Couldn't find %s" % path)
            print("json file paths should be relative to paths_cfg")
            sys.exit(1)

        json_paths[direction] = path

    if not figpath:
        figpath = os.path.join(indir, "nbr_grid.%s" % format)
        log_file_path = os.path.join(indir, "nbr_grid.log")
    else:
        figpath = util.abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    LOGGER.log_file_path = log_file_path
    plot_data = {}
    for direction, path in list(json_paths.items()):
        LOGGER.input_file(path)
        data = util.load_loglin_stats(path)
        plot_data[direction] = data

    fig = draw_position_grid(plot_data, sample_size)

    fig.text(0.4, 0.955, "Ending Base", fontsize=20)
    fig.text(0.03, 0.55, "Starting Base", rotation=90, fontsize=20)
    fig.tight_layout(rect=(0.06, 0, 0.95, 0.95))

    fig.savefig(figpath)
    LOGGER.output_file(figpath)
    print("Wrote", figpath)


def draw_spectrum_grid(data, plot_cfg=None, sample_size=False,
                       width=8, height=8,
                       title_space=1.1,
                       axis_font_size=20, tick_font_size=10, ylim=None):
    all_bases = 'ACGT'
    bases = list(data)
    bases.sort()
    assert set(data.keys()) <= set(all_bases)
    num_bases = len(data)
    figsize = list(plot_cfg.get('grid', 'figsize'))
    if num_bases == 2:
        figsize[1] = figsize[1] // 2

    f, arr = pyplot.subplots(num_bases, 4, sharex=True,
                             sharey=True, figsize=figsize)
    f.subplots_adjust(hspace=.2, wspace=0.2)
    if ylim is None:
        ylim = 0
        for base in bases:
            ylim = max(ylim, max(data[base].values()))

        ylim *= 1.2

    major_loc = MultipleLocator(ylim / 2.)
    yformat = FuncFormatter(format_float(1e-3, float_places=2))

    for i in range(num_bases):
        a = arr[i][0]
        r = a.yaxis.set_major_locator(major_loc)
        r = a.yaxis.set_major_formatter(yformat)
        a.set_ylim([0, ylim])

    xtick_labels = []
    xlim = None
    x, y = None, 0
    dx = None
    letter_pad = 0.1
    for i in range(num_bases):
        for j in range(4):
            a = arr[i][j]
            r = [a.tick_params(axis='x', which='major',
                               bottom='off', top='off')]
            xtick_labels.append(a.get_xticklabels())
            xlim = xlim or a.get_xlim()
            if x is None and xlim:
                x = letter_pad * xlim[1]
                dx = xlim[1] - (2 * letter_pad)

            if i == j:
                continue

            dy = data[bases[i]][all_bases[j]]
            if dy < 0:
                invert = True
                dy *= -1
            else:
                invert = False

            text.add_letter(all_bases[j], x, y, dx, dy, ax=a, invert=invert)

        pyplot.setp(xtick_labels, visible=False)

    for i in range(num_bases):
        a = arr[i][0]
        a.set_ylabel(bases[i], rotation=0, labelpad=20, fontsize=24)

    for i in range(4):
        a = arr[0][i]
        a.set_title(all_bases[i], fontsize=24)

    if num_bases == 4:
        rect = (0.06, 0, 0.95, 0.95)
        title_y = 0.955
    elif num_bases == 2:
        rect = (0.06, 0, 0.95, 0.9)
        title_y = 0.9
    else:
        raise ValueError("Number of bases must be 4 or 2, not %d" % num_bases)

    f.text(0.4, title_y, "Ending Base", fontsize=20)
    f.text(0.03, 0.55, "Starting Base", rotation=90, fontsize=20)
    f.tight_layout(rect=rect)

    return f


def load_spectra_data(json_path, group_col):
    # for each starting base, we need the total relative entropy
    # we need the ret's for each ending base
    LOGGER.input_file(json_path)
    data = util.load_loglin_stats(json_path)
    bases = list(data)
    bases.sort()
    assert set(data.keys()) <= set('CTAG')

    if group_col:
        assert group_col in "strand group", \
            "group_col must be 'group' or 'strand', got %s" % group_col

    if 'group' in data[bases[0]]['stats'].columns:
        group_col = 'group'
    else:
        group_col = 'strand'

    selected_group = {'strand': '+', 'group': '1'}.get(group_col, None)
    assert selected_group is not None, selected_group

    result = {}
    for base in bases:
        total_re = data[base]['rel_entropy']
        subset = data[base]['stats'][data[base]['stats'][
            group_col].apply(str) == selected_group].copy()
        if subset.empty:
            print("No entries equal to '%s'" % str(selected_group))
            exit(-1)

        total_ret = numpy.fabs(subset["ret"]).sum()
        subset['prop'] = total_re * subset['ret'] / total_ret
        subset['end'] = [d[-1:] for d in subset['direction']]
        result[base] = dict((b, v)
                            for i, b, v in
                            subset[['end', 'prop']].to_records())

    return result


_json_path = click.option('--json_path', required=True,
              help="Path to spectra analysis spectra_analysis.json")
_group_label = click.option('--group_label', default='group', help="Id for reference group")


@main.command()
@_json_path
@_group_label
@_plot_cfg
@_no_type3
@_figpath
@_format
@_sample_size
@_force_overwrite
@_dry_run
def spectra_grid(json_path, group_label, plot_cfg, no_type3, figpath, format, sample_size,
                 force_overwrite, dry_run):
    """draws logo from mutation spectra analysis"""
    # the following is for logging
    args = locals()
    if no_type3:
        util.exclude_type3_fonts()

    if not figpath:
        dirname = os.path.dirname(json_path)
        figpath = os.path.join(dirname, "spectra_grid.%s" % format)
        log_file_path = os.path.join(dirname, "spectra_grid.log")
    else:
        figpath = util.abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    LOGGER.log_file_path = log_file_path

    LOGGER.log_message(str(args), label='vars')

    data = load_spectra_data(json_path, group_label)

    if plot_cfg:
        LOGGER.input_file(plot_cfg)
    plot_cfg = util.get_plot_configs(cfg_path=plot_cfg)
    f = draw_spectrum_grid(data, sample_size=sample_size, plot_cfg=plot_cfg)
    f.savefig(figpath)
    LOGGER.output_file(figpath)
    print("Wrote", figpath)


@main.command()
@_json_path
@_plot_cfg
@_no_type3
@_figpath
@_format
@_sample_size
@_force_overwrite
@_dry_run
def mi(json_path, plot_cfg, no_type3, figpath, format, sample_size,
       force_overwrite, dry_run):
    """draws conventional sequence logo, using MI, from first order effects"""
    # the following is for logging
    json_path = util.abspath(json_path)
    args = locals()
    if no_type3:
        util.exclude_type3_fonts()

    if not figpath:
        dirname = os.path.dirname(json_path)
        figpath = os.path.join(dirname, "MI.%s" % format)
        log_file_path = os.path.join(dirname, "MI.log")
    else:
        figpath = util.abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    LOGGER.log_file_path = log_file_path

    if plot_cfg:
        LOGGER.input_file(plot_cfg)

    LOGGER.log_message(str(args), label='vars')

    data = util.load_loglin_stats(json_path)
    positions = list(data.keys())
    positions.sort()
    num_pos = len(positions) + 1
    mp = num_pos // 2
    counts_array = numpy.zeros((4, num_pos), int)
    for i, pos in enumerate(positions):
        if i >= mp:
            i += 1
        pos_stats = data[pos]['stats']
        counts = pos_stats[pos_stats['mut'] == 'M'][["base", "count"]]
        counts = dict(zip(counts['base'], counts['count']))
        for base in counts:
            base_index = DNA.alphabet.index(base)
            counts_array[base_index, i] = counts[base]

    freq_matrix = entropy.counts_to_freq_matrix(counts_array)
    mit = entropy.get_mit(freq_matrix, freq_matrix=True)
    mi = mit.sum(axis=0)
    char_hts = get_mi_char_heights(numpy.fabs(mit), mi)

    plot_cfg = util.get_plot_configs(cfg_path=plot_cfg)
    figsize = plot_cfg.get('1-way plot', 'figsize')
    ytick_font = plot_cfg.get('1-way plot', 'ytick_fontsize')
    xtick_font = plot_cfg.get('1-way plot', 'xtick_fontsize')
    ylabel_font = plot_cfg.get('1-way plot', 'ylabel_fontsize')
    xlabel_font = plot_cfg.get('1-way plot', 'xlabel_fontsize')
    fig = logo.draw_multi_position(char_hts.T,
                                   characters=[list(DNA)] * num_pos,
                                   position_indices=list(range(num_pos)),
                                   figsize=figsize,
                                   figwidth=figsize[0],
                                   xtick_fontsize=xtick_font,
                                   ytick_fontsize=ytick_font,
                                   sort_data=True)

    ax = fig.gca()
    ax.tick_params(axis='y', labelsize=ytick_font)
    ax.tick_params(axis='x', labelsize=xtick_font)
    ax.set_ylabel("MI", fontsize=ylabel_font)
    ax.set_xlabel("Position", fontsize=xlabel_font)
    fig.tight_layout()
    fig.savefig(figpath)
    LOGGER.output_file(figpath)
    print("Wrote", figpath)
