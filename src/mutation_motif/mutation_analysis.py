#!/usr/bin/env python
import os
from itertools import combinations

import click
from scitrack import CachingLogger

from cogent3 import make_table
from cogent3.maths.stats import chisqprob
from mutation_motif import draw, log_lin, motif_count, spectra_analysis, util

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-2020, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


LOGGER = CachingLogger(create_dir=True)


def make_summary(results):
    """returns records from analyses as list"""
    rows = []
    for position_set in results:
        if type(position_set) != str:
            position = ":".join(position_set)
        else:
            position = position_set

        re = results[position_set]["rel_entropy"]
        dev = results[position_set]["deviance"]
        df = results[position_set]["df"]
        prob = results[position_set]["prob"]
        formula = results[position_set]["formula"]
        rows.append([position, re, dev, df, prob, formula])
    return rows


def get_grouped_combined_counts(table, position, group_label):
    """wraps motif_count.get_combined_counts for groups"""
    group_cats = table.distinct_values(group_label)
    all_data = []
    header = None
    for category in group_cats:
        subtable = table.filtered(lambda x: x == category, columns=group_label)
        counts = motif_count.get_combined_counts(subtable, position)
        if header is None:
            header = [group_label] + list(counts.header)

        counts = counts.with_new_column(
            group_label, lambda x: category, columns=counts.header[0]
        )
        all_data.extend(counts.tolist(header))
    counts = make_table(header=header, rows=all_data)
    counts.sorted(columns=[group_label, "mut"])
    return counts


def get_position_effects(table, position_sets, group_label=None):
    pos_results = {}
    grouped = group_label is not None
    if grouped:
        assert len(table.distinct_values(group_label)) == 2

    for position_set in position_sets:
        if not grouped:
            counts = motif_count.get_combined_counts(table, position_set)
        else:
            counts = get_grouped_combined_counts(
                table, position_set, group_label=group_label
            )
        rel_entropy, deviance, df, stats, formula = log_lin.position_effect(
            counts, group_label=group_label
        )
        if deviance < 0:
            p = 1.0
        else:
            p = chisqprob(deviance, df)

        pos_results[position_set] = dict(
            rel_entropy=rel_entropy,
            deviance=deviance,
            df=df,
            stats=stats,
            formula=formula,
            prob=p,
        )
    return pos_results


def single_position_effects(table, positions, group_label=None):
    single_results = get_position_effects(table, positions, group_label=group_label)
    return single_results


def get_two_position_effects(table, positions, group_label=None):
    two_pos_results = get_position_effects(
        table, list(combinations(positions, 2)), group_label=group_label
    )
    return two_pos_results


def get_three_position_effects(table, positions, group_label=None):
    three_pos_results = get_position_effects(
        table, list(combinations(positions, 3)), group_label=group_label
    )
    return three_pos_results


def get_four_position_effects(table, positions, group_label=None):
    result = get_position_effects(
        table, list(combinations(positions, 4)), group_label=group_label
    )
    return result


def single_group(
    counts_table, outpath, group_label, group_ref, positions, first_order, dry_run,
):
    # Collect statistical analysis results
    summary = []

    max_results = {}
    # Single position analysis
    print("Doing single position analysis")
    single_results = single_position_effects(
        counts_table, positions, group_label=group_label
    )
    summary += make_summary(single_results)

    max_results[1] = max(single_results[p]["rel_entropy"] for p in single_results)
    if not dry_run:
        outfilename = os.path.join(outpath, "1.json")
        util.dump_loglin_stats(single_results, outfilename)
        LOGGER.output_file(outfilename, label="analysis1")

    fig = draw.get_1way_position_drawable(
        single_results, None, group_label=group_label, group_ref=group_ref
    )

    if not dry_run:
        outfilename = os.path.join(outpath, "1.pdf")
        fig.write(outfilename)
        LOGGER.output_file(outfilename)

    if first_order:
        msg = "Done! Check %s for your results" % outpath
        summary = make_table(
            header=["Position", "RE", "Deviance", "df", "prob", "formula"],
            rows=summary,
            digits=2,
            space=2,
        )
        if not dry_run:
            outfilename = os.path.join(outpath, "summary.txt")
            summary.write(outfilename, sep="\t")
            LOGGER.output_file(outfilename, label="summary")

        return msg

    print("Doing two positions analysis")
    results = get_two_position_effects(counts_table, positions, group_label=group_label)
    summary += make_summary(results)

    max_results[2] = max(results[p]["rel_entropy"] for p in results)
    if not dry_run:
        outfilename = os.path.join(outpath, "2.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis2")

    fig = draw.get_2way_position_drawable(
        results, None, group_label=group_label, group_ref=group_ref
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "2.pdf")
        fig.write(outfilename)
        LOGGER.output_file(outfilename)

    print("Doing three positions analysis")
    results = get_three_position_effects(
        counts_table, positions, group_label=group_label
    )
    summary += make_summary(results)

    max_results[3] = max(results[p]["rel_entropy"] for p in results)
    if not dry_run:
        outfilename = os.path.join(outpath, "3.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis3")

    fig = draw.get_3way_position_drawable(
        results, None, group_label=group_label, group_ref=group_ref,
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "3.pdf")
        fig.write(outfilename)
        LOGGER.output_file(outfilename)

    print("Doing four positions analysis")
    results = get_four_position_effects(
        counts_table, positions, group_label=group_label
    )
    summary += make_summary(results)

    max_results[4] = max(results[p]["rel_entropy"] for p in results)
    if not dry_run:
        outfilename = os.path.join(outpath, "4.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis4")

    fig = draw.get_4way_position_drawable(
        results, None, group_label=group_label, group_ref=group_ref,
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "4.pdf")
        fig.write(outfilename)
        LOGGER.output_file(outfilename)

    # now generate summary plot
    summary = make_table(
        header=["Position", "RE", "Deviance", "df", "prob", "formula"],
        rows=summary,
        digits=2,
        space=2,
    )
    fig = draw.get_summary_drawable(summary, None)
    if not dry_run:
        outfilename = os.path.join(outpath, "summary.pdf")
        fig.write(outfilename)
        LOGGER.output_file(outfilename)

    if not dry_run:
        outfilename = os.path.join(outpath, "summary.txt")
        summary.write(outfilename, sep="\t")
        LOGGER.output_file(outfilename, label="summary")

    msg = "Done! Check %s for your results" % outpath
    return msg


_countsfile = click.option("-1", "--countsfile", help="tab delimited file of counts.")
_outpath = click.option("-o", "--outpath", help="Directory path to write data.")
_countsfile2 = click.option(
    "-2", "--countsfile2", help="second group motif counts file."
)
_strand_symmetry = click.option(
    "-s",
    "--strand_symmetry",
    is_flag=True,
    help="single counts file but second group is strand.",
)
_force_overwrite = click.option(
    "-F", "--force_overwrite", is_flag=True, help="Overwrite existing files."
)
_dry_run = click.option(
    "-D",
    "--dry_run",
    is_flag=True,
    help="Do a dry run of the analysis without writing output.",
)
_verbose = click.option("-v", "--verbose", is_flag=True, help="Display more output.")


@click.group()
def main():
    pass


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


@main.command()
@_countsfile
@_outpath
@_countsfile2
@_first_order
@_strand_symmetry
@_group_label
@_group_ref
@_verbose
@_dry_run
def nbr(
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
            group_label, lambda x: "1", columns=counts_table.header[0]
        )

        fn2 = util.abspath(countsfile2)
        counts_table2 = util.load_table_from_delimited_file(fn2, sep="\t")

        LOGGER.input_file(fn2, label="countsfile2_path")

        counts_table2 = counts_table2.with_new_column(
            group_label, lambda x: "2", columns=counts_table2.header[0]
        )
        # now combine
        header = [group_label] + counts_table2.header[:-1]
        raw1 = counts_table1.tolist(header)
        raw2 = counts_table2.tolist(header)
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
        counts_table, outpath, group_label, group_ref, positions, first_order, dry_run,
    )
    click.secho(msg, fg="green")


@main.command()
@_countsfile
@_outpath
@_countsfile2
@_strand_symmetry
@_force_overwrite
@_dry_run
@_verbose
def spectra(
    countsfile,
    outpath,
    countsfile2,
    strand_symmetry,
    force_overwrite,
    dry_run,
    verbose,
):
    """log-linear analysis of mutation spectra between groups
    """
    spectra_analysis.main(
        countsfile,
        outpath,
        countsfile2,
        strand_symmetry,
        force_overwrite,
        dry_run,
        verbose,
    )
    click.secho("Done spectra!", fg="green")
