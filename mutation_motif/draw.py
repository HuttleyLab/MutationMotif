import os, sys
from optparse import make_option
from ConfigParser import SafeConfigParser

import numpy
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter

from cogent import LoadTable, DNA
from cogent.util.option_parsing import parse_command_line_parameters

from mutation_motif import util, mutation_analysis, logo
from mutation_motif.height import get_re_char_heights

from scitrack import CachingLogger

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

def get_selected_indices(stats): # from mutation_analysis
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
        mut_stats = mut_stats.sort('ret')
        characters[index] = list(mut_stats['base'])
        rets[:,index] = mut_stats['ret']

    heights = get_re_char_heights(rets, re_positionwise=position_re)
    return heights.T, characters, range(num_pos)

def draw_position_grid(directions, sample_size=False, width=8, height=8, title_space=1.1, axis_font_size=20, tick_font_size=10, ylim=None):
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
            positions = data.keys()
            positions.sort()
        number = data[positions[0]]['stats']["count"].sum() / 2
        heights, characters, indices = get_plot_data(data, positions)
        adaptive_y = max(adaptive_y, logo.est_ylim(heights))
        plottables.append([direction, heights, characters, indices, number])
    
    if ylim is None:
        ylim = adaptive_y
    
    for direction, heights, characters, indices, number in plottables:
        fr, to = map(bases.index, direction.split('to'))
        ax = axes[fr, to]
        fig = logo.draw_multi_position(heights, characters=characters,
                        position_indices=indices, ylim=ylim, ax=ax, figwidth=width, verbose=False)
        if sample_size:
            y = ax.get_ylim()[1]
            ax.text(0.2, y * 0.85, "N={:,}".format(number), fontsize=10)
    
    xformat = FuncFormatter(format_float(1e-3, float_places=2))

    for i in range(4):
        top_ax = axes[0, i]
        top_ax.set_title(bases[i], fontsize=axis_font_size, weight="bold", y=1.1)
        
        lft_ax = axes[i, 0]
        for yticklabel in lft_ax.get_yticklabels():
            yticklabel.set_fontsize(tick_font_size)
            yticklabel.set_rotation(0)
        lft_ax.yaxis.set_major_formatter(FuncFormatter(xformat))
        lft_ax.set_ylabel(bases[i], rotation=0, fontsize=axis_font_size, weight="bold")
        lft_ax.yaxis.labelpad = axis_font_size
        
        btm_ax = axes[-1, i]
        for xticklabel in btm_ax.get_xticklabels():
            xticklabel.set_fontsize(tick_font_size)
            xticklabel.set_rotation(0)
    
    f.tight_layout()
    
    return f

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "\n".join([
"draw grid of neighbour effect logos",
"",
"Takes the 1.json files produced by mutation_nbr for all mutation directions."])

script_info['required_options'] = [
     make_option('-p','--paths_cfg',
                 help='Text file listing path for 1.json file for each '\
                 'mutation direction (e.g. AtoG).'),
    ]

script_info['optional_options'] = [
     make_option('--figpath', help='Filename for grid plot. Overides format.'),
    make_option('--format', default='pdf', choices=['pdf', 'png'],
        help='Plot format [default:%default].'),
    make_option('--sample_size', action='store_true', default=False,
        help='Include sample size on each subplot.'),
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
    
    config_path = util.abspath(opts.paths_cfg)
    indir = os.path.dirname(config_path)
    parser = SafeConfigParser()
    parser.optionxform = str # stops automatic conversion to lower case
    parser.read(config_path)
    
    json_paths = {}
    for direction, path in parser.items("json_paths"):
        # assumes paths are relative to indir
        path = os.path.join(indir, path)
        if not os.path.exists(path):
            print "Couldn't find %s" % path
            print "json file paths should be relative to paths_cfg"
            sys.exit(1)
        
        json_paths[direction] = path
    
    if not opts.figpath:
        figpath = os.path.join(indir, "mutation_grid.%s" % opts.format)
        log_file_path = os.path.join(indir, "mutation_grid_draw.log")
    else:
        figpath = util.abspath(opts.figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])
    
    LOGGER.log_file_path = log_file_path
    plot_data = {}
    for direction, path in json_paths.items():
        LOGGER.input_file(path)
        data = util.load_loglin_stats(path)
        plot_data[direction] = data
    
    fig = draw_position_grid(plot_data, opts.sample_size)
    
    fig.savefig(figpath)
    LOGGER.output_file(figpath)
    print "Wrote", figpath

if __name__ == "__main__":
    main()
    