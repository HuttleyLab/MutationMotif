import glob
import os
import sys
import time
from collections import OrderedDict
from collections.abc import Mapping
from configparser import ConfigParser
from importlib import resources

import click
from cogent3 import make_table
from scitrack import CachingLogger

from mutation_motif import __version__, spectra_analysis, util
from mutation_motif.all_counts import check_found_filenames, direction
from mutation_motif.aln_to_counts import align_to_counts, get_counts_filename
from mutation_motif.complement import make_strand_symmetric_table
from mutation_motif.draw import (
    get_1way_position_drawable,
    get_2way_position_drawable,
    get_3way_position_drawable,
    get_4way_position_drawable,
    get_grid_drawable,
    get_mi_plot_data,
    get_position_grid_drawable,
    get_spectra_grid_drawable,
    get_summary_drawable,
    write_pdf,
)
from mutation_motif.mutation_analysis import single_group
from mutation_motif.util import (
    abspath,
    get_nbr_path_config,
    get_subtables,
    load_loglin_stats,
    load_table_from_delimited_file,
    makedirs,
)

_click_command_opts = dict(
    no_args_is_help=True,
    context_settings={"show_default": True},
)
_countsfile = click.option("-1", "--countsfile", help="tab delimited file of counts.")
_outpath = click.option("-o", "--outpath", help="Directory path to write data.")
_countsfile2 = click.option(
    "-2",
    "--countsfile2",
    help="second group motif counts file.",
)
_strand_symmetry = click.option(
    "-s",
    "--strand_symmetry",
    is_flag=True,
    help="single counts file but second group is strand.",
)
_force_overwrite = click.option(
    "-F",
    "--force_overwrite",
    is_flag=True,
    help="Overwrite existing files.",
)
_dry_run = click.option(
    "-D",
    "--dry_run",
    is_flag=True,
    help="Do a dry run of the analysis without writing output.",
)
_verbose = click.option("-v", "--verbose", is_flag=True, help="Display more output.")
_first_order = click.option(
    "--first_order",
    is_flag=True,
    help="Consider only first order effects. Defaults "
    "to considering up to 4th order interactions.",
)
_group_label = click.option("-g", "--group_label", help="second group label.")
_group_ref = click.option(
    "-r",
    "--group_ref",
    default=None,
    help="reference group value for results presentation.",
)
_paths_cfg = click.option(
    "--paths_cfg",
    required=True,
    help="Text file listing path for 1.json file for "
    "each mutation direction (e.g. AtoG).",
)
_figpath = click.option(
    "--figpath",
    help="Filename for plot file. Suffix defines format.",
)
_plot_cfg = click.option(
    "--plot_cfg",
    help="Config file for plot size, font size settings.",
)
_fig_cfg = click.option("--fig_config", required=True, type=click.Path(exists=True))


class OrderedGroup(click.Group):
    """custom group class to ensure help function returns commands in desired order.
    class is adapted from Максим Стукало's answer to
    https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
    """

    def __init__(
        self,
        name: str | None = None,
        commands: Mapping[str, click.Command] | None = None,
        **kwargs,
    ):
        super().__init__(name, commands, **kwargs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx: click.Context) -> Mapping[str, click.Command]:
        return self.commands


@click.group(cls=OrderedGroup)
@click.version_option(__version__)  # add version option
def main():
    """mm: point mutation analysis tools"""


@main.command(**_click_command_opts)
@click.option(
    "-a",
    "--align_path",
    required=True,
    help="fasta aligned file centred on mutated position.",
)
@click.option("-o", "--output_path", required=True, help="Path to write data.")
@click.option(
    "-f",
    "--flank_size",
    required=True,
    type=int,
    help="Number of bases per side to include.",
)
@click.option(
    "--direction",
    required=True,
    default=None,
    type=click.Choice(
        [
            "AtoC",
            "AtoG",
            "AtoT",
            "CtoA",
            "CtoG",
            "CtoT",
            "GtoA",
            "GtoC",
            "GtoT",
            "TtoA",
            "TtoC",
            "TtoG",
        ],
    ),
    help="Mutation direction.",
)
@click.option(
    "-S",
    "--seed",
    help="Seed for random number generator (e.g. 17, or 2015-02-13)."
    " Defaults to system time.",
)
@click.option(
    "-R",
    "--randomise",
    is_flag=True,
    help="Randomises the observed data, observed and reference "
    "counts distributions should match.",
)
@click.option(
    "--step",
    default="1",
    type=click.Choice(["1", "2", "3"]),
    help='Specifies a "frame" for selecting the random base.',
)
@_dry_run
@_force_overwrite
def prep_nbr(
    align_path,
    output_path,
    flank_size,
    direction,
    seed,
    randomise,
    step,
    dry_run,
    force_overwrite,
):
    """export tab delimited counts table from alignment centred on SNP position

    Output file is written to the same path with just the file suffix changed
    from fasta to txt."""
    LOGGER = CachingLogger(create_dir=True)
    if not seed:
        seed = str(time.time())

    LOGGER.log_args()

    align_path = abspath(align_path)
    output_path = abspath(output_path)

    counts_filename = get_counts_filename(align_path, output_path)
    runlog_path = counts_filename.replace(".txt", ".log")
    LOGGER.log_file_path = runlog_path

    if not dry_run:
        if not force_overwrite and (
            os.path.exists(counts_filename) or os.path.exists(runlog_path)
        ):
            msg = (
                "Either %s or %s already exist. Force overwrite of existing"
                " files with -F."
            )
            raise ValueError(msg % (counts_filename, runlog_path))

        makedirs(output_path)

        LOGGER.input_file(align_path, label="align_path")
        LOGGER.log_message(str(seed), label="random_number_seed")

    start_time = time.time()

    # run the program

    counts_table = align_to_counts(
        align_path,
        output_path,
        flank_size,
        direction,
        step,
        seed,
        randomise,
        dry_run,
        LOGGER=LOGGER,
    )

    if not dry_run:
        counts_table.write(counts_filename, sep="\t")
        LOGGER.output_file(counts_filename)

    # determine runtime
    duration = time.time() - start_time
    if not dry_run:
        LOGGER.log_message(f"{duration / 60.0:.2f}", label="run duration (minutes)")
    LOGGER.shutdown()


@main.command(**_click_command_opts)
@click.option(
    "-c",
    "--counts_pattern",
    help="glob pattern uniquely identifying all 12 mutation counts files.",
)
@click.option("-o", "--output_path", help="Path to write combined_counts data.")
@click.option(
    "-s",
    "--strand_symmetric",
    is_flag=True,
    help="produces table suitable for strand symmetry test.",
)
@click.option(
    "-p",
    "--split_dir",
    default=None,
    help="path to write individual direction strand symmetric tables.",
)
@_dry_run
@_force_overwrite
def prep_spectra(
    counts_pattern,
    output_path,
    strand_symmetric,
    split_dir,
    dry_run,
    force_overwrite,
):
    """export tab delimited combined counts table by appending the 12 mutation
    direction tables, adding a new column ``direction``."""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    output_path = abspath(output_path)
    split_dir = abspath(split_dir) if strand_symmetric and split_dir else None
    # check we the glob pattern produces the correct number of files
    counts_files = glob.glob(counts_pattern)
    check_found_filenames(counts_files)

    counts_filename = os.path.join(output_path, "combined_counts.txt")
    runlog_path = os.path.join(output_path, "combined_counts.log")

    if not dry_run:
        if not force_overwrite and (
            os.path.exists(counts_filename) or os.path.exists(runlog_path)
        ):
            msg = (
                "Either %s or %s already exist. Force overwrite of "
                "existing files with -F."
            )
            raise ValueError(msg % (counts_filename, runlog_path))

        makedirs(output_path)
        if split_dir:
            makedirs(split_dir)

        LOGGER.log_file_path = runlog_path
        for fn in counts_files:
            LOGGER.input_file(fn, label="count_file")

    start_time = time.time()

    # run the program
    basenames = []
    tables = []
    for fn in counts_files:
        basenames.append(os.path.basename(fn))
        mutation = direction.findall(fn)[0]
        table = load_table_from_delimited_file(fn, sep="\t")
        table.title = mutation
        tables.append(table)

    table = tables[0].appended("direction", tables[1:])

    if strand_symmetric:
        table = make_strand_symmetric_table(table)

    if split_dir:
        group_subtables = get_subtables(table, group_label="direction")

    if not dry_run:
        table.write(counts_filename, sep="\t")
        LOGGER.output_file(counts_filename)

        if split_dir:
            for group, subtable in group_subtables:
                # we first assume that group is part of the filenames!
                fn = [bn for bn in basenames if group in bn]
                fn = fn[0] if len(fn) == 1 else f"{group}.txt"
                counts_filename = os.path.join(split_dir, fn)
                subtable.write(counts_filename, sep="\t")
                LOGGER.output_file(counts_filename)

    # determine runtime
    duration = time.time() - start_time
    if not dry_run:
        LOGGER.log_message(f"{duration / 60.0:.2f}", label="run duration (minutes)")

    click.secho("Done!", fg="green")
    LOGGER.shutdown()


@main.command(**_click_command_opts)
@_countsfile
@_outpath
@_countsfile2
@_first_order
@_strand_symmetry
@_group_label
@_group_ref
@_verbose
@_dry_run
def ll_nbr(
    countsfile,
    outpath,
    countsfile2,
    first_order,
    strand_symmetry,
    group_label,
    group_ref,
    verbose,
    dry_run,
):
    """log-linear analysis of neighbouring base influence on point mutation

    Writes estimated statistics, figures and a run log to the specified
    directory outpath.

    See documentation for count table format requirements.
    """
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    outpath = util.abspath(outpath)

    if not dry_run:
        util.makedirs(outpath)
        runlog_path = os.path.join(outpath, "analysis.log")
        LOGGER.log_file_path = runlog_path

    counts_filename = util.abspath(countsfile)
    counts_table = util.load_table_from_delimited_file(counts_filename, sep="\t")

    LOGGER.input_file(counts_filename, label="countsfile1_path")

    positions = [c for c in counts_table.header if c.startswith("pos")]
    if not first_order and len(positions) != 4:
        raise ValueError("Requires four positions for analysis")

    group_label = group_label or None
    group_ref = group_ref or None
    if strand_symmetry:
        group_label = "strand"
        group_ref = group_ref or "+"
        if group_label not in counts_table.header:
            print("ERROR: no column named 'strand', exiting.")
            exit(-1)

    if countsfile2:
        print("Performing 2 group analysis")
        group_label = group_label or "group"
        group_ref = group_ref or "1"
        counts_table1 = counts_table.with_new_column(
            group_label,
            lambda x: "1",
            columns=counts_table.header[0],
        )

        fn2 = util.abspath(countsfile2)
        counts_table2 = util.load_table_from_delimited_file(fn2, sep="\t")

        LOGGER.input_file(fn2, label="countsfile2_path")

        counts_table2 = counts_table2.with_new_column(
            group_label,
            lambda x: "2",
            columns=counts_table2.header[0],
        )
        # now combine
        header = [group_label] + list(counts_table2.header[:-1])
        raw1 = counts_table1.to_list(header)
        raw2 = counts_table2.to_list(header)
        counts_table = make_table(header=header, rows=raw1 + raw2)

        if not dry_run:
            outfile = os.path.join(outpath, "group_counts_table.txt")
            counts_table.write(outfile, sep="\t")
            LOGGER.output_file(outfile, label="group_counts")

    if dry_run or verbose:
        print()
        print(counts_table)
        print()

    msg = single_group(
        counts_table,
        outpath,
        group_label,
        group_ref,
        positions,
        first_order,
        dry_run,
        LOGGER,
    )
    click.secho(msg, fg="green")


@main.command(**_click_command_opts)
@_countsfile
@_outpath
@_countsfile2
@_strand_symmetry
@_force_overwrite
@_dry_run
@_verbose
def ll_spectra(
    countsfile,
    outpath,
    countsfile2,
    strand_symmetry,
    force_overwrite,
    dry_run,
    verbose,
):
    """log-linear analysis of mutation spectra between groups"""
    LOGGER = CachingLogger(create_dir=True)
    spectra_analysis.main(
        countsfile,
        outpath,
        countsfile2,
        strand_symmetry,
        force_overwrite,
        dry_run,
        verbose,
        LOGGER,
    )
    click.secho("Done spectra!", fg="green")


@main.command()
@_paths_cfg
@_plot_cfg
@_figpath
@_force_overwrite
@_dry_run
def draw_nbr_matrix(
    paths_cfg,
    plot_cfg,
    figpath,
    force_overwrite,
    dry_run,
):
    """draws square matrix of sequence logo's from neighbour analysis"""
    LOGGER = CachingLogger(create_dir=True)
    config_path = abspath(paths_cfg)
    indir = os.path.dirname(config_path)
    if not figpath:
        figpath = os.path.join(indir, "nbr_matrix.%s" % format)
        log_file_path = os.path.join(indir, "nbr_matrix.log")
    else:
        figpath = abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    if not force_overwrite and os.path.exists(figpath):
        click.secho(f"{figpath} alreadyt exists")
        sys.exit(0)

    LOGGER.log_args()

    parser = ConfigParser()
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

    LOGGER.log_file_path = log_file_path
    plot_data = {}
    for direction, path in json_paths.items():
        LOGGER.input_file(path)
        data = load_loglin_stats(path)
        plot_data[direction] = path

    fig = get_position_grid_drawable(plot_data, plot_cfg)

    write_pdf(fig, figpath)
    LOGGER.output_file(figpath)
    click.secho(f"Wrote {figpath}", fg="green")
    LOGGER.shutdown()


@main.command()
@_fig_cfg
@_figpath
def draw_grid(fig_config, figpath):
    """draws an arbitrary shaped grid of mutation motifs based on a config file"""
    # we read in the config file and determine number of rows and columns
    # paths, headings, etc ..
    # then create the figure and axes and call the mutation_motif drawing code
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    if not figpath:
        dirname = os.path.dirname(fig_config.name)
        figpath = os.path.join(dirname, "drawn_array.%s" % format)
        log_file_path = os.path.join(dirname, "drawn_array.log")
    else:
        figpath = abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    makedirs(os.path.dirname(figpath))
    LOGGER.log_file_path = log_file_path

    fig = get_grid_drawable(fig_config)
    write_pdf(fig, figpath)
    click.secho(f"Wrote {figpath}", fg="green")
    LOGGER.shutdown()


_json_path = click.option(
    "-p",
    "--json_path",
    required=True,
    help="Path to spectra analysis spectra_analysis.json",
)
_fig_cfg = click.option("--fig_config", required=True, type=click.Path(exists=True))


@main.command()
@_json_path
@_group_label
@_plot_cfg
@_figpath
@_force_overwrite
@_dry_run
def draw_spectra_grid(
    json_path,
    group_label,
    plot_cfg,
    figpath,
    force_overwrite,
    dry_run,
):
    """draws logo from mutation spectra analysis"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    if not figpath:
        dirname = os.path.dirname(json_path)
        figpath = os.path.join(dirname, "spectra_grid.%s" % format)
        log_file_path = os.path.join(dirname, "spectra_grid.log")
    else:
        figpath = abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    LOGGER.log_file_path = log_file_path

    # data = load_spectra_data(json_path, group_label)

    if plot_cfg:
        LOGGER.input_file(plot_cfg)

    fig = get_spectra_grid_drawable(
        json_path,
        plot_cfg=plot_cfg,
        group_label=group_label,
    )
    write_pdf(fig, figpath)
    LOGGER.output_file(figpath)
    click.secho(f"Wrote {figpath}", fg="green")
    LOGGER.shutdown()


@main.command()
@click.option(
    "-p",
    "--json_paths",
    type=click.Path(exists=True),
    help="config file with json paths",
)
@_plot_cfg
@_group_label
@_force_overwrite
@_dry_run
def draw_nbr(
    json_paths,
    plot_cfg,
    group_label,
    force_overwrite,
    dry_run,
):
    """makes motifs for independent or higher order interactions"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    dirname = os.path.dirname(json_paths)
    LOGGER.log_file_path = os.path.join(dirname, "nbr.log")

    if plot_cfg:
        LOGGER.input_file(plot_cfg)

    paths = get_nbr_path_config(json_paths)
    one_way = "1-way plot"
    two_way = "2-way plot"
    three_way = "3-way plot"
    four_way = "4-way plot"
    funcs = {
        one_way: get_1way_position_drawable,
        two_way: get_2way_position_drawable,
        three_way: get_3way_position_drawable,
        four_way: get_4way_position_drawable,
    }

    for order in (one_way, two_way, three_way, four_way):
        if order not in paths:
            continue
        LOGGER.input_file(paths[order].inpath)
        data = load_loglin_stats(paths[order].inpath)
        fig = funcs[order](data, plot_cfg, group_label=group_label)
        write_pdf(fig, paths[order].outpath)
        LOGGER.output_file(paths[order].outpath)
        click.secho(f"Wrote {paths[order].outpath}", fg="green")

    summary = "summary"
    if summary in paths:
        LOGGER.input_file(paths[summary].inpath)
        fig = get_summary_drawable(paths[summary].inpath, plot_cfg)
        write_pdf(fig, paths[summary].outpath)
        LOGGER.output_file(paths[summary].outpath)
        click.secho(f"Wrote {paths[summary].outpath}", fg="green")

    click.secho("Done!", fg="green")
    LOGGER.shutdown()


@main.command()
@_json_path
@_plot_cfg
@_group_label
@_figpath
@click.option(
    "--use_freq",
    is_flag=True,
    help="Use freqs rather than residuals for letter height.",
)
@_force_overwrite
@_dry_run
def draw_mi(
    json_path,
    plot_cfg,
    group_label,
    figpath,
    use_freq,
    force_overwrite,
    dry_run,
):
    """draws conventional sequence logo, using MI, from first order effects"""
    LOGGER = CachingLogger(create_dir=True)
    global mi_use_freqs
    mi_use_freqs = use_freq

    LOGGER.log_args()
    # the following is for logging
    json_path = abspath(json_path)

    if not figpath:
        dirname = os.path.dirname(json_path)
        figpath = os.path.join(dirname, "MI.pdf")
        log_file_path = os.path.join(dirname, "MI.log")
    else:
        figpath = abspath(figpath)
        log_file_path = "%s.log" % ".".join(figpath.split(".")[:-1])

    LOGGER.log_file_path = log_file_path

    if plot_cfg:
        LOGGER.input_file(plot_cfg)

    data = load_loglin_stats(json_path)
    fig = get_1way_position_drawable(
        data,
        plot_cfg,
        group_label=group_label,
        get_heights=get_mi_plot_data,
    )
    for ann in fig.layout.annotations:
        if "RE" in ann.text:
            ann.text = ann.text.replace("RE", "MI")
            break
    write_pdf(fig, figpath)
    LOGGER.output_file(figpath)
    click.secho(f"Wrote {figpath}", fg="green")
    LOGGER.shutdown()


@main.command()
@click.argument("outpath")
def draw_export_cfg(outpath):
    """exports the sample config files to the nominated path"""
    import shutil

    if os.path.exists(outpath):
        click.secho(
            "outpath already exists, delete it or choose different dest",
            fg="red",
        )
        sys.exit(1)

    path = resources.files("mutation_motif") / "cfgs"
    shutil.copytree(path, outpath)
    click.secho(f"Contents written to {outpath}", fg="green")


_sample_size = click.option(
    "--sample_size",
    is_flag=True,
    help="Include sample size on each subplot.",
)
