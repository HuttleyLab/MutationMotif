import os
from itertools import combinations

from cogent3 import make_table

from mutation_motif import (
    draw,
    log_lin,
    motif_count,
    util,
)


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
            group_label,
            lambda x: category,
            columns=counts.header[0],
        )
        all_data.extend(counts.to_list(header))
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
                table,
                position_set,
                group_label=group_label,
            )
        result = log_lin.position_effect(
            counts,
            group_label=group_label,
        )
        p = result.pvalue
        pos_results[position_set] = dict(
            rel_entropy=result.relative_entropy,
            deviance=result.deviance,
            df=result.nfp,
            stats=result.df,
            formula=result.formula,
            prob=p,
        )
    return pos_results


def single_position_effects(table, positions, group_label=None):
    return get_position_effects(table, positions, group_label=group_label)


def get_two_position_effects(table, positions, group_label=None):
    return get_position_effects(
        table,
        list(combinations(positions, 2)),
        group_label=group_label,
    )


def get_three_position_effects(table, positions, group_label=None):
    return get_position_effects(
        table,
        list(combinations(positions, 3)),
        group_label=group_label,
    )


def get_four_position_effects(table, positions, group_label=None):
    return get_position_effects(
        table,
        list(combinations(positions, 4)),
        group_label=group_label,
    )


def single_group(
    counts_table,
    outpath,
    group_label,
    group_ref,
    positions,
    first_order,
    dry_run,
    LOGGER,
):
    # Collect statistical analysis results
    summary = []

    max_results = {}
    # Single position analysis
    print("Doing single position analysis")
    single_results = single_position_effects(
        counts_table,
        positions,
        group_label=group_label,
    )
    summary += make_summary(single_results)

    max_results[1] = max(single_results[p]["rel_entropy"] for p in single_results)
    if not dry_run:
        outfilename = os.path.join(outpath, "1.json")
        util.dump_loglin_stats(single_results, outfilename)
        LOGGER.output_file(outfilename, label="analysis1")

    fig = draw.get_1way_position_drawable(
        single_results,
        None,
        group_label=group_label,
        group_ref=group_ref,
    )
    write_fig = util.pdf_writer()
    if not dry_run:
        outfilename = os.path.join(outpath, "1.pdf")
        write_fig(fig, outfilename)
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
        results,
        None,
        group_label=group_label,
        group_ref=group_ref,
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "2.pdf")
        write_fig(fig, outfilename)
        LOGGER.output_file(outfilename)

    print("Doing three positions analysis")
    results = get_three_position_effects(
        counts_table,
        positions,
        group_label=group_label,
    )
    summary += make_summary(results)

    max_results[3] = max(results[p]["rel_entropy"] for p in results)
    if not dry_run:
        outfilename = os.path.join(outpath, "3.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis3")

    fig = draw.get_3way_position_drawable(
        results,
        None,
        group_label=group_label,
        group_ref=group_ref,
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "3.pdf")
        write_fig(fig, outfilename)
        LOGGER.output_file(outfilename)

    print("Doing four positions analysis")
    results = get_four_position_effects(
        counts_table,
        positions,
        group_label=group_label,
    )
    summary += make_summary(results)

    max_results[4] = max(results[p]["rel_entropy"] for p in results)
    if not dry_run:
        outfilename = os.path.join(outpath, "4.json")
        util.dump_loglin_stats(results, outfilename)
        LOGGER.output_file(outfilename, label="analysis4")

    fig = draw.get_4way_position_drawable(
        results,
        None,
        group_label=group_label,
        group_ref=group_ref,
    )
    if not dry_run:
        outfilename = os.path.join(outpath, "4.pdf")
        write_fig(fig, outfilename)
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
        write_fig(fig, outfilename)
        LOGGER.output_file(outfilename)

    if not dry_run:
        outfilename = os.path.join(outpath, "summary.txt")
        summary.write(outfilename, sep="\t")
        LOGGER.output_file(outfilename, label="summary")

    msg = "Done! Check %s for your results" % outpath
    LOGGER.shutdown()
    return msg
