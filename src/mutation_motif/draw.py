import os
import sys
from configparser import ConfigParser

import click
import numpy
from pkg_resources import resource_filename
from scitrack import CachingLogger

from cogent3.core.profile import MotifCountsArray
from cogent3.draw.drawable import Drawable, get_domain
from cogent3.draw.logo import get_base_logo_layout, get_logo
from cogent3.draw.logo import get_mi_char_heights as c3_get_mi
from cogent3.util.union_dict import UnionDict
from mutation_motif.height import get_mi_char_heights, get_re_char_heights
from mutation_motif.util import (abspath, est_ylim, get_grid_config,
                                 get_nbr_config, get_nbr_matrix_config,
                                 get_nbr_path_config,
                                 get_order_max_re_from_summary,
                                 get_position_number, get_selected_indices,
                                 get_spectra_config, get_summary_config,
                                 load_loglin_stats, makedirs)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-2020, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "BSD-3"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


LOGGER = CachingLogger(create_dir=True)


_axis_lines = dict(
    mirror=True,
    linewidth=1,
    showgrid=False,
    linecolor="black",
    showline=True,
    visible=True,
    zeroline=False,
)
COLOURS = dict(A="green", C="blue", G="orange", T="red")

mi_use_freqs = False


def get_mi_plot_data(pwise_results, positions, group_label=None, group_ref=None):
    """calculate character heights from counts data based on mutual information

    Parameters
    ----------
    pwise_results : dict
        statistics from log-linear modelling using the mutation_analysis
        functions
    positions :  list/tuple
        series of positions names, e.g. ['pos0', 'pos1'].
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.

    Returns
    -------
    [[(base1, height1), (base2, height2), ..], ...]. The middle position has
    an empty dict.
    """
    mid = (len(positions) + 1) // 2

    base_order = list("CTAG")
    counts = []
    for pos in positions:
        stats = pwise_results[pos]["stats"]
        mut_stats = stats[
            get_selected_indices(stats, group_label=group_label, group_ref=group_ref)
        ][["base", "count"]]
        arr = mut_stats.to_numpy()
        c = arr[:, 1].astype(int)
        b = arr[:, 0]
        data = dict(zip(b, c))
        counts.append([data[b] for b in base_order])
    counts = MotifCountsArray(counts, base_order)
    freqs = counts.to_freq_array()
    if not mi_use_freqs:
        mit = 2 - freqs.entropy_terms().array
        mi = 2 - freqs.entropy()

        char_heights = []
        for i in range(mi.shape[0]):
            hts = get_mi_char_heights(mit[i], mi[i])
            char_heights.append({b: v for b, v in zip(base_order, hts)})
    else:
        char_heights = []
        hts = c3_get_mi(freqs)
        for i in range(freqs.shape[0]):
            char_heights.append(hts[i].to_dict())

    char_heights.insert(mid, {})
    return char_heights


def get_re_plot_data(pwise_results, positions, group_label=None, group_ref=None):
    """calculate character heights based on relative entropy

    Parameters
    ----------
    pwise_results : dict
        statistics from log-linear modelling using the mutation_analysis
        functions
    positions :  list/tuple
        series of positions names, e.g. ['pos0', 'pos1'].
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.

    Returns
    -------
    [[(base1, height1), (base2, height2), ..], ...]. The middle position has
    an empty dict.

    Notes
    -----

    The statistics are derived from the log-linear modelling fit.
    """
    mid = (len(positions) + 1) // 2

    char_heights = []
    for index, pos in enumerate(positions):
        stats = pwise_results[pos]["stats"]
        mut_stats = stats[
            get_selected_indices(stats, group_label=group_label, group_ref=group_ref)
        ][["base", "ret"]]
        arr = mut_stats.to_numpy()
        rets_ = arr[:, 1].astype(float)
        chars_ = arr[:, 0]
        hts = get_re_char_heights(
            rets_, re_positionwise=pwise_results[pos]["rel_entropy"]
        )

        vals = {b: ret for b, ret in zip(chars_, hts)}
        char_heights.append(vals)

    char_heights.insert(mid, {})
    return char_heights


def _get_nbr_matrix_axnum(bases):
    """returns {(row, col): axis number, ...}"""
    coords_to_axis = {}
    dir_coords = {}
    n = 1
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            coords_to_axis[(i, j)] = n
            n += 1
            dir_coords[f"{bases[i]}to{bases[j]}"] = (i, j)

    return coords_to_axis, dir_coords


def get_matrix_row_col_titles(
    labels,
    total=None,
    coord=None,
    is_row=True,
    textangle=0,
    font_size=14,
    pad=0,
    space=0.05,
):
    """construct plotly text annotations for row/column labels

    Parameters
    ----------

    labels : series
        elements are strings that will be displayed
    total : int or None
        total number of blocks on axis. If provided, the number of labels must be
        <= this. Can not provide coord too.
    coord : float or None
        coord is the y-axis coordinate if is_row, otherwise it's the x-axis coord
    is_row : bool
        If True, positioning is on the y-axis.
    textangle : float
        Orientation of the text. Defaults to 0.
    font_size : int

    pad : int
        In pixels, defines the shift relative to coordinate along the axis
        defined by is_row.
    space : float
        Spacing between rows/columns of plot. Used in conjunction with total to
        compute text placement in the middle of a row / column.

    Returns
    -------

    [{annot1}, {annot2}, ...]
    """
    assert not (total and coord), "provide one of"
    if total:
        assert len(labels) <= total, "Number of labels must be <= total"

    titles = []
    for index in range(len(labels)):
        label = labels[index]
        title = UnionDict(
            {
                "font": {"size": font_size},
                "showarrow": False,
                "text": f"<b>{label}</b>",
                "xref": "paper",
                "yref": "paper",
                "xanchor": "center",
                "yanchor": "middle",
                "textangle": textangle,
            }
        )
        # note that plotly display is in cartesian, so the y array coordinate
        # needs to be reversed
        if total:
            coord = sum(get_domain(total, index, is_y=is_row, space=space)) / 2
        if is_row:
            title.x = 0
            title.y = coord
            title.xshift = pad
        else:
            title.x = coord
            title.y = 1
            title.yanchor = "bottom"
            title.yshift = pad
        titles.append(title)
    return titles


def get_position_grid_drawable(directions, plot_cfg, ylim=None):
    """produces 4x4 Drawable of mutation motifs for all mutation direction for
    1st order effects.

    Parameters
    ----------
    directions : dict
        Keys are the mutation directions in 'XtoY' form, where X, Y are nucleotides.
        Values can be either a path to a json file containing results from a
        log-linear analysis, or the pre-loaded contents of such a file.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    default_cfg = get_nbr_matrix_config(None)
    if plot_cfg:
        # over-ride with user settings
        plot_cfg = get_grid_config(plot_cfg)
        default_cfg |= plot_cfg

    plot_cfg = default_cfg

    # needs to use a plot config
    bases = list("CTAG")
    layout = UnionDict(shapes=[], annotations=[])
    coords_to_axis, dir_to_coords = _get_nbr_matrix_axnum(bases)
    sep = plot_cfg.get("space", 0.05)
    positions = None
    plottables = UnionDict()
    for direction in directions:
        value = directions[direction]
        if not isinstance(value, dict):
            data = load_loglin_stats(value)
        else:
            data = value

        if positions is None:
            positions = list(data.keys())
            positions.sort()
        heights = get_re_plot_data(data, positions)
        plottables[direction] = heights

    if ylim is None and "ylim" not in plot_cfg:
        all_heights = []
        for data in plottables.values():
            for h in data:
                all_heights.append(sum(abs(v) for v in h.values()))
        ylim = max(0, est_ylim(all_heights))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    # customised xticks, text and location
    mid = (len(positions) + 1) // 2
    xtick_vals = [i for i in range(len(positions) + 1)]
    xtick_text = [f"{i - mid}" for i in xtick_vals]

    for direction in plottables:
        row, col = dir_to_coords[direction]
        axnum = coords_to_axis[(row, col)]
        # setup subplot layout
        base_layout = get_base_logo_layout(
            axnum, plot_cfg.xtick_fontsize, plot_cfg.ytick_fontsize
        )
        layout |= base_layout
        ax = "" if axnum == 1 else f"{axnum}"
        # domains
        layout[f"xaxis{ax}"]["domain"] = get_domain(4, col, is_y=False, space=sep)
        # note that plotly display is cartesian, so the y array coordinate
        # needs to be reversed
        layout[f"yaxis{ax}"]["domain"] = get_domain(4, row, is_y=True, space=sep)
        # tweaks to axis display
        layout[f"xaxis{ax}"].tickvals = xtick_vals
        layout[f"xaxis{ax}"].ticktext = xtick_text
        layout[f"yaxis{ax}"].title = None
        if col != 0 and (row, col) != (0, 1):
            layout[f"yaxis{ax}"] |= dict(showticklabels=False)

        if row != 3 and (row, col) != (2, 3):
            layout[f"xaxis{ax}"] |= dict(showticklabels=False)

        r = get_logo(plottables[direction], ylim=ylim, axnum=axnum, layout=base_layout)
        # cogent3 returns a Drawable instance, so we access it's layout
        r = r.figure.layout
        layout.shapes.extend(r.shapes)

    # configure row, column titles and turn off ticks
    row_titles = get_matrix_row_col_titles(
        bases,
        4,
        is_row=True,
        font_size=plot_cfg.ylabel_fontsize,
        pad=plot_cfg.ylabel_pad,
    )
    col_titles = get_matrix_row_col_titles(
        bases,
        4,
        is_row=False,
        font_size=plot_cfg.xlabel_fontsize,
        pad=plot_cfg.xlabel_pad,
    )
    layout.annotations.extend(col_titles + row_titles)
    # outer title Ending Base & Starting Base
    col_titles = get_matrix_row_col_titles(
        ["Ending Base"],
        1,
        is_row=False,
        pad=plot_cfg.col_title_pad,
        font_size=plot_cfg.col_title_fontsize,
    )
    row_titles = get_matrix_row_col_titles(
        ["Starting Base"],
        1,
        textangle=270,
        is_row=True,
        pad=plot_cfg.row_title_pad,
        font_size=plot_cfg.row_title_fontsize,
    )
    layout.annotations.extend(row_titles + col_titles)

    layout.template = "plotly_white"
    layout.margin = plot_cfg.margin
    fig = Drawable(layout=layout, width=plot_cfg.width, height=plot_cfg.height)

    return fig


# get effect order summary stuff
def get_summary_drawable(data, plot_cfg, ylim=None):
    """produces bar chart Drawable of maximum RE for all effect orders

    Parameters
    ----------
    data : str, Table
        path to a tsv file generated by the nbr command, or a cogent3 Table
        produced from reading such a file.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    default_cfg = get_summary_config(None)
    if plot_cfg:
        plot_cfg = get_summary_config(plot_cfg)
        default_cfg |= plot_cfg

    plot_cfg = default_cfg

    order_vals = get_order_max_re_from_summary(data)
    order, stat = list(zip(*order_vals))

    if ylim is None and "ylim" not in plot_cfg:
        ylim = (0, est_ylim(stat))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    trace = UnionDict(y=stat, x=order, type="bar")
    layout = UnionDict(template="plotly_white")
    layout.annotations = get_matrix_row_col_titles(
        ["Effect Order"],
        coord=0.5,
        is_row=False,
        pad=plot_cfg.xlabel_pad,
        font_size=plot_cfg.xlabel_fontsize,
    )
    layout.annotations[0]["y"] = 0
    layout.annotations += get_matrix_row_col_titles(
        ["RE<sub><i>max</i></sub>"],
        coord=0.5,
        is_row=True,
        textangle=270,
        pad=plot_cfg.ylabel_pad,
        font_size=plot_cfg.ylabel_fontsize,
    )
    layout.yaxis = UnionDict(
        range=[0, ylim],
        tickfont=dict(size=plot_cfg.ytick_fontsize),
        titlefont=dict(size=plot_cfg.ylabel_fontsize),
    )
    layout.xaxis = UnionDict(
        tickmode="array",
        ticktext=list([str(o) for o in order]),
        tickvals=order,
        tickfont=dict(size=plot_cfg.xtick_fontsize),
        titlefont=dict(size=plot_cfg.xlabel_fontsize),
    )

    axis_lines = dict(
        mirror=True,
        linewidth=1,
        showgrid=False,
        linecolor="black",
        showline=True,
        visible=True,
        zeroline=False,
    )
    layout.xaxis |= axis_lines
    layout.yaxis |= axis_lines
    layout.margin = plot_cfg.margin
    height = plot_cfg.height
    width = plot_cfg.width
    fig = Drawable(layout=layout, height=height, width=width, showlegend=False)
    fig.add_trace(trace)
    return fig


# get independent position drawable, a single plot
def get_1way_position_drawable(
    data,
    plot_cfg,
    group_label=None,
    group_ref=None,
    ylim=None,
    get_heights=get_re_plot_data,
):
    """produces mutation motif Drawable for a single mutation direction of
    independent position effects.

    Parameters
    ----------
    data : dict
        Statistics from loading a 1.json file produced by the nbr command.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.
    get_heights : callable
        Function reference used to compute character heights from log-linear
        analysis statistics. Function must satisfy signature defined by
        get_re_plot_data.

    Returns
    -------
    cogent3 Drawable.
    """
    plot_cfg = get_nbr_config(plot_cfg, "1-way plot")
    positions = sorted(
        [k for k in data if k.startswith("pos")], key=get_position_number
    )
    heights = get_heights(data, positions, group_label=group_label, group_ref=group_ref)
    if ylim is None and "ylim" not in plot_cfg:
        values = []
        for h in heights:
            values.append(sum(abs(v) for v in h.values()))
        ylim = max(0, est_ylim(values))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    fig = get_logo(heights, ylim=ylim)
    layout = fig.figure.layout
    layout.yaxis.title = None
    # customised xticks, text and location
    mid = (len(positions) + 1) // 2
    xtick_vals = [i for i in range(len(positions) + 1)]
    xtick_text = [f"{i - mid}" for i in xtick_vals]
    layout.xaxis.tickvals = xtick_vals
    layout.xaxis.ticktext = xtick_text

    layout.annotations = get_matrix_row_col_titles(
        ["Position"],
        coord=0.5,
        is_row=False,
        pad=plot_cfg.xlabel_pad,
        font_size=plot_cfg.xlabel_fontsize,
    )
    layout.annotations[0]["y"] = 0
    layout.annotations += get_matrix_row_col_titles(
        ["RE"],
        coord=0.5,
        is_row=True,
        textangle=270,
        pad=plot_cfg.ylabel_pad,
        font_size=plot_cfg.ylabel_fontsize,
    )
    layout.margin = plot_cfg.margin
    layout.template = "plotly_white"
    height = plot_cfg.height
    width = plot_cfg.width
    fig = Drawable(layout=layout, height=height, width=width)
    return fig


# get multi-way interactions drawable
def _get_multi_way_position_drawables(
    data,
    plot_cfg,
    group_coords,
    group_domains,
    bases,
    group_label=None,
    group_ref=None,
    ylim=None,
):
    """produces mutation motif Drawable for orders >= 2.

    Parameters
    ----------
    data : dict
        Statistics from loading a <n>.json file produced by the nbr command where
        n is either 2, 3 or 4. Keys are tuples of position names.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    group_coords : dict
        Key's match those of the data dict. Values correspond to x, y placement.
    group_domains : dict
        Key's match those of the data dict. Values correspond to [start, end] of
        the plotly figure domain.
    bases : series
        column labels corresponding to log-linear analysis result tables of the
        positions
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.
    get_heights : callable
        Function reference used to compute character heights from log-linear
        analysis statistics. Function must satisfy signature defined by
        get_re_plot_data.

    Returns
    -------
    cogent3 Drawable.
    """
    bases = list(bases)
    layout = UnionDict(annotations=[], shapes=[])

    if ylim is None and "ylim" not in plot_cfg:
        heights = [data[group]["rel_entropy"] for group in data]
        ylim = max(0, est_ylim(heights))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    positions = {
        m: get_position_number(m)
        for group in data
        for m in group
        if m.startswith("pos")
    }
    num_pos = len(positions) + 1
    mid_pos = num_pos // 2
    ordered = [""] * num_pos
    for label, index in positions.items():
        if index >= mid_pos:
            index += 1
        ordered[index] = label

    # customised xticks, text and location
    mid = (len(positions) + 1) // 2
    xtick_vals = [i for i in range(len(positions) + 1)]
    xtick_text = [f"{i - mid}" for i in xtick_vals]

    axnum = 1

    num_row = num_col = 0
    for group in sorted(data):
        col, row = group_coords[group]
        num_row = max(num_row, row)
        num_col = max(num_col, col)

    for group in sorted(data):
        char_heights = [None] * num_pos
        col, row = group_coords[group]

        base_layout = get_base_logo_layout(
            axnum, plot_cfg.xtick_fontsize, plot_cfg.ytick_fontsize
        )
        layout |= base_layout
        ax = "" if axnum == 1 else f"{axnum}"
        # domains
        layout[f"xaxis{ax}"] |= dict(
            domain=group_domains[group].x,
            tickvals=xtick_vals,
            ticktext=xtick_text,
            range=[-0.5, num_pos + 0.25],
        )
        # note that plotly display is cartesian, so the y array coordinate
        # needs to be reversed
        layout[f"yaxis{ax}"] |= dict(
            domain=group_domains[group].y, title=None, range=[0, ylim]
        )
        # make sure xaxis has correct range

        if col != 0:
            layout[f"yaxis{ax}"] |= dict(showticklabels=False)

        if row != num_row:
            layout[f"xaxis{ax}"] |= dict(showticklabels=False)

        indices = list(map(ordered.index, group))
        rel_entropy = data[group]["rel_entropy"]
        stats = data[group]["stats"]
        idx = get_selected_indices(stats, group_label=group_label, group_ref=group_ref)
        mut_stats = stats[idx][bases + ["ret"]]
        hts = get_re_char_heights(list(mut_stats["ret"]), re_positionwise=rel_entropy)

        for i, base in enumerate(bases):
            chars = list(mut_stats[base])
            char_re = list(zip(chars, hts))
            char_heights[indices[i]] = char_re

        r = get_logo(char_heights, axnum=axnum, ylim=ylim, layout=base_layout)
        r = r.figure.layout
        layout.shapes.extend(r.shapes)

        axnum += 1

    # add x/y axis labels as annotations
    xlabel = get_matrix_row_col_titles(
        ["Position"],
        coord=0.5,
        is_row=False,
        pad=plot_cfg.xlabel_pad,
        font_size=plot_cfg.xlabel_fontsize,
    )
    xlabel[0].y = 0
    ylabel = get_matrix_row_col_titles(
        ["RE"],
        coord=0.5,
        is_row=True,
        textangle=270,
        pad=plot_cfg.ylabel_pad,
        font_size=plot_cfg.ylabel_fontsize,
    )
    layout.annotations = xlabel + ylabel
    height = plot_cfg.height
    width = plot_cfg.width
    layout.margin = plot_cfg.margin
    layout.template = "plotly_white"
    fig = Drawable(layout=layout, height=height, width=width)
    return fig


# get 2-way interactions drawable, a lower triangular grid plot
def get_2way_position_drawable(
    data, plot_cfg, group_label=None, group_ref=None, ylim=None
):
    """produces mutation motif Drawable for a single mutation direction of
    2nd order interactions.

    Parameters
    ----------
    data : dict
        Statistics from loading a 2.json file produced by the nbr command.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    plot_cfg = get_nbr_config(plot_cfg, "2-way plot")
    sep = plot_cfg.get("space", 0.03)
    positions = {
        m: get_position_number(m)
        for group in data
        for m in group
        if m.startswith("pos")
    }
    group_coords = {(a, b): (positions[a], positions[b]) for a, b in sorted(data)}
    domains = UnionDict()
    for pair, (col, row) in group_coords.items():
        y = get_domain(len(positions) - 1, row - 1, is_y=True, space=sep)
        x = get_domain(len(positions) - 1, col, is_y=False, space=sep)
        domains[pair] = UnionDict(x=x, y=y)

    fig = _get_multi_way_position_drawables(
        data,
        plot_cfg,
        group_coords,
        domains,
        ("base1", "base2"),
        group_label=group_label,
        group_ref=group_ref,
        ylim=ylim,
    )
    return fig


# get 3-way interactions drawable, a lower triangular grid plot
def get_3way_position_drawable(
    data, plot_cfg, group_label=None, group_ref=None, ylim=None
):
    """produces mutation motif Drawable for a single mutation direction of
    3rd order interactions.

    Parameters
    ----------
    data : dict
        Statistics from loading a 2.json file produced by the nbr command.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    plot_cfg = get_nbr_config(plot_cfg, "3-way plot")
    sep = plot_cfg.get("space", 0.03)
    num_rows = plot_cfg.num_rows
    group_coords = {
        triple: divmod(index, num_rows) for index, triple in enumerate(sorted(data))
    }
    domains = UnionDict()
    for group, (col, row) in group_coords.items():
        y = get_domain(num_rows, row, is_y=True, space=sep)
        x = get_domain(num_rows, col, is_y=False, space=sep)
        domains[group] = UnionDict(x=x, y=y)

    fig = _get_multi_way_position_drawables(
        data,
        plot_cfg,
        group_coords,
        domains,
        ("base1", "base2", "base3"),
        group_label=group_label,
        group_ref=group_ref,
        ylim=ylim,
    )
    return fig


def get_4way_position_drawable(
    data, plot_cfg, group_label=None, group_ref=None, ylim=None
):
    """produces mutation motif Drawable for a single mutation direction of
    4th order interactions.

    Parameters
    ----------
    data : dict
        Statistics from loading a 2.json file produced by the nbr command.
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    plot_cfg = get_nbr_config(plot_cfg, "4-way plot")
    num_rows = plot_cfg.get("num_rows", 1)
    group_coords = {
        four: divmod(index, num_rows) for index, four in enumerate(sorted(data))
    }
    domains = UnionDict()
    for group, (col, row) in group_coords.items():
        y = get_domain(num_rows, row, is_y=True, space=0.03)
        x = get_domain(num_rows, col, is_y=False, space=0.03)
        domains[group] = UnionDict(x=x, y=y)

    fig = _get_multi_way_position_drawables(
        data,
        plot_cfg,
        group_coords,
        domains,
        ("base1", "base2", "base3", "base4"),
        group_label=group_label,
        group_ref=group_ref,
        ylim=ylim,
    )
    return fig


def get_spectra_row(data, bases, axnum, ylim=1, domain=(0, 1), colours=None):
    """computes single row of mutation spectra plot

    Parameters
    ----------
    data : dict
        {base: relative entropy term, }
    bases
        ordered list of bases
    axnum : int
        the axis number, unique to this row
    ylim : float
        limit of y-axis
    domain
        list of start, end of vertical plotm domain spanned by this row
    colours : dict
        maps bases to fill colors

    Returns
    -------
    UnionDict in plotly layout structure
    """
    from cogent3.draw.letter import get_character

    axis_lines = _axis_lines.copy()

    layout = UnionDict(shapes=[])
    xref = "x" if axnum == 1 else f"x{axnum}"
    yref = "y" if axnum == 1 else f"y{axnum}"
    anchor = "" if axnum == 1 else axnum
    layout["xaxis"] = dict(
        domain=[0, 1], showticklabels=False, range=[1, 5], anchor=f"y{anchor}"
    )
    layout["yaxis"] = dict(
        domain=domain, range=[0, ylim], anchor=f"x{anchor}", ticks="inside"
    )
    layout["xaxis"] |= axis_lines
    layout["yaxis"] |= axis_lines

    for i, b in enumerate(bases):
        if b not in data:
            continue

        val = data[b]
        c = get_character(
            letter=b, x=i + 1, y=0, height=val, width=1, fillcolor=colours[b]
        )
        c.xref = xref
        c.yref = yref
        layout.shapes.append(c.as_shape())

    return layout


def get_spectra_grid_drawable(
    data,
    plot_cfg=None,
    bases=("A", "C", "G", "T"),
    ylim=None,
    group_label=None,
    group_ref=None,
):
    """computes mutation spectra Drawable from log-linear analysis

    Parameters
    ----------
    data : dict
        nested dicts of from/to bases and ending relative entropy term
    bases
        base order
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable
    """
    # get default cfg values
    default_cfg = get_spectra_config(None)
    if plot_cfg:
        plot_cfg = get_spectra_config(plot_cfg)
        default_cfg |= plot_cfg

    plot_cfg = default_cfg

    if group_ref is None:
        group_ref = {"strand": "+", "group": "1"}.get(group_label, None)

    assert group_ref is not None, group_ref

    data = load_spectra_data(data, group_label=group_label, group_ref=group_ref)

    if ylim is None and "ylim" not in plot_cfg:
        heights = [abs(data[a][b]) for a in data for b in data[a]]
        ylim = max(0, est_ylim(heights))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    sep = plot_cfg.get("space", 0.04)
    assert set(data.keys()) <= set(bases)
    # arrange the subset of letters in bases order
    ordered = sorted(data, key=bases.index)

    colours = COLOURS

    # a subplot for each top level character
    layout = UnionDict(shapes=[], annotations=[], template="plotly_white")

    for i, b in enumerate(ordered):
        axnum = i + 1
        axis_name = "axis" if axnum == 1 else f"axis{axnum}"
        domain = get_domain(len(ordered), i, is_y=True, space=sep)
        r = get_spectra_row(
            data[b],
            bases=bases,
            axnum=axnum,
            ylim=ylim,
            domain=domain,
            colours=colours,
        )

        layout[f"x{axis_name}"] = r.xaxis
        layout[f"y{axis_name}"] = r.yaxis
        layout.shapes.extend(r.shapes)

    # get column and row titles
    col_titles = get_matrix_row_col_titles(
        bases,
        total=4,
        is_row=False,
        pad=plot_cfg.xlabel_pad,
        font_size=plot_cfg.xlabel_fontsize,
    )
    row_titles = get_matrix_row_col_titles(
        ordered,
        total=4,
        is_row=True,
        pad=plot_cfg.ylabel_pad,
        font_size=plot_cfg.ylabel_fontsize,
    )
    layout.annotations.extend(col_titles + row_titles)

    # outer title Ending Base & Starting Base
    col_titles = get_matrix_row_col_titles(
        [plot_cfg.col_title],
        coord=0.5,
        is_row=False,
        pad=plot_cfg.col_title_pad,
        font_size=plot_cfg.col_title_fontsize,
    )
    row_titles = get_matrix_row_col_titles(
        [plot_cfg.row_title],
        coord=0.5,
        textangle=270,
        is_row=True,
        pad=plot_cfg.row_title_pad,
        font_size=plot_cfg.row_title_fontsize,
    )
    layout.annotations.extend(col_titles + row_titles)
    layout.margin = plot_cfg.margin
    fig = Drawable(layout=layout, width=plot_cfg.width, height=plot_cfg.height)
    return fig


def load_spectra_data(json_path, group_label, group_ref):
    """load results from log-linear mutation spectrum analysis

    Parameters
    ----------
    json_path : str
        path to file produced by spectra command
    group_label : str
        Factor name, e.g. "group", from which counts are extracted.
    group_ref : str
        Factor value, e.g. "1", from which counts are extracted.

    Returns
    -------
    dict of {(from, to): relative entropy, ...}
    """
    # for each starting base, we need the total relative entropy
    # we need the ret's for each ending base
    if isinstance(json_path, dict):
        data = json_path
    else:
        data = load_loglin_stats(json_path)

    bases = list(data)
    bases.sort()
    assert set(data.keys()) <= set("CTAG")

    if group_label:
        assert group_label in "strand group", (
            "group_label must be 'group' or 'strand', got %s" % group_label
        )

    if "group" in data[bases[0]]["stats"].columns:
        group_label = "group"
    else:
        group_label = "strand"

    result = {}
    for base in bases:
        total_re = data[base]["rel_entropy"]
        subset = data[base]["stats"][
            data[base]["stats"][group_label].apply(str) == group_ref
        ].copy()
        if subset.empty:
            print("No entries equal to '%s'" % str(group_ref))
            exit(-1)

        total_ret = numpy.fabs(subset["ret"]).sum()
        subset["prop"] = total_re * subset["ret"] / total_ret
        subset["end"] = [d[-1:] for d in subset["direction"]]
        result[base] = dict((b, v) for i, b, v in subset[["end", "prop"]].to_records())

    return result


_json_path = click.option(
    "-p",
    "--json_path",
    required=True,
    help="Path to spectra analysis " "spectra_analysis.json",
)
_group_label = click.option("--group_label", help="Id for reference group")


_fig_cfg = click.option("--fig_config", required=True, type=click.Path(exists=True))


def _get_axis_num(nrow, row, col):
    """
    computes axis number

    Parameters
    ----------
    nrow : int
        number of rows
    row : int
        row number
    col : int
        column number
    """
    return (nrow * row + col) + 1


def get_grid_drawable(plot_cfg, ylim=None):
    """produces mutation motif Drawable for a 1st order position effects for
    arbitrary shaped grid.

    Parameters
    ----------
    plot_cfg : str, None
        Path to a config file, or None. In latter case, uses built-in defaults.
        Use the export_cfg command to get sample cfg files..
    ylim : float
        Maximum value for all y-axis plots.

    Returns
    -------
    cogent3 Drawable.
    """
    # get default cfg values
    default_cfg = get_grid_config(None)
    if plot_cfg:
        # over-ride with user settings
        plot_cfg = get_grid_config(plot_cfg)
        default_cfg |= plot_cfg
        plot_cfg = default_cfg

    layout = UnionDict(shapes=[], annotations=[])
    plottable = {}
    coords_to_axis = {}
    positions = None
    for coord in plot_cfg.paths:
        data = load_loglin_stats(plot_cfg.paths[coord])
        coords_to_axis[tuple(coord)] = _get_axis_num(
            plot_cfg.num_rows, coord[1], coord[0]
        )
        if positions is None:
            positions = list(data)
            positions.sort()
        else:
            assert set(positions) == set(data), "positions must match"
        heights = get_re_plot_data(data, positions)
        plottable[coord] = heights

    if ylim is None and "ylim" not in plot_cfg:
        all_heights = []
        for data in plottable.values():
            for h in data:
                all_heights.append(sum(abs(v) for v in h.values()))
        ylim = max(0, est_ylim(all_heights))
    elif "ylim" in plot_cfg:
        ylim = plot_cfg.ylim

    # customised xticks, text and location
    mid = (len(positions) + 1) // 2
    xtick_vals = [i for i in range(len(positions) + 1)]
    xtick_text = [f"{i - mid}" for i in xtick_vals]

    for coord in plottable:
        col, row = coord
        axnum = coords_to_axis[tuple(coord)]
        base_layout = get_base_logo_layout(
            axnum, plot_cfg.xtick_fontsize, plot_cfg.ytick_fontsize
        )
        layout |= base_layout
        ax = "" if axnum == 1 else f"{axnum}"
        # domains
        layout[f"xaxis{ax}"]["domain"] = get_domain(
            plot_cfg.num_cols, col, is_y=False, space=plot_cfg.space
        )
        # note that plotly display is cartesian, so the y array coordinate
        # needs to be reversed
        layout[f"yaxis{ax}"]["domain"] = get_domain(
            plot_cfg.num_rows, row, is_y=True, space=plot_cfg.space
        )
        # tweaks to axis display
        layout[f"xaxis{ax}"].tickvals = xtick_vals
        layout[f"xaxis{ax}"].ticktext = xtick_text
        layout[f"yaxis{ax}"].title = None
        if col != 0:
            layout[f"yaxis{ax}"] |= dict(showticklabels=False)

        if row != plot_cfg.num_rows - 1:
            layout[f"xaxis{ax}"] |= dict(showticklabels=False)

        r = get_logo(plottable[coord], axnum=axnum, ylim=ylim, layout=base_layout)
        # cogent3 returns a Drawable instance, so we access it's layout
        r = r.figure.layout
        layout.shapes.extend(r.shapes)

    if "row_titles" in plot_cfg:
        row_titles = get_matrix_row_col_titles(
            plot_cfg.row_titles,
            plot_cfg.num_rows,
            is_row=True,
            font_size=plot_cfg.row_title_fontsize,
            pad=plot_cfg.row_title_pad,
        )
        layout.annotations.extend(row_titles)
    if "col_titles" in plot_cfg:
        col_titles = get_matrix_row_col_titles(
            plot_cfg.col_titles,
            plot_cfg.num_cols,
            is_row=False,
            font_size=plot_cfg.col_title_fontsize,
            pad=plot_cfg.col_title_pad,
        )
        layout.annotations.extend(col_titles)

    x_title = {
        "font": {"size": plot_cfg.xlabel_fontsize},
        "showarrow": False,
        "text": "Position",
        "x": 0.5,
        "xanchor": "center",
        "xref": "paper",
        "y": 0,
        "yshift": plot_cfg.xlabel_pad,
        "yanchor": "bottom",
        "yref": "paper",
    }

    y_title = {
        "font": {"size": plot_cfg.ylabel_fontsize},
        "showarrow": False,
        "text": "RE",
        "x": 0,
        "xshift": plot_cfg.ylabel_pad,
        "xanchor": "center",
        "xref": "paper",
        "y": 0.5,
        "yanchor": "bottom",
        "yref": "paper",
    }
    layout.annotations.extend([x_title, y_title])
    layout.template = "plotly_white"
    layout.margin = plot_cfg.margin

    fig = Drawable(layout=layout, width=plot_cfg.width, height=plot_cfg.height)
    return fig


# the CLI functions and options
@click.group()
def main():
    """draw mutation motif logo's and spectra"""
    pass


# defining the CLI options
_paths_cfg = click.option(
    "--paths_cfg",
    required=True,
    help="Text file listing path for 1.json file for "
    "each mutation direction (e.g. AtoG).",
)

_figpath = click.option(
    "--figpath", help="Filename for plot file. Suffix defines format."
)
_plot_cfg = click.option(
    "--plot_cfg", help="Config file for plot size, font size settings."
)

_sample_size = click.option(
    "--sample_size", is_flag=True, help="Include sample size on each subplot."
)

_force_overwrite = click.option(
    "-F", "--force_overwrite", is_flag=True, help="Overwrite existing files."
)
_dry_run = click.option(
    "-D",
    "--dry_run",
    is_flag=True,
    help="Do a dry run of the analysis without writing " "output.",
)


@main.command()
@_paths_cfg
@_plot_cfg
@_figpath
@_force_overwrite
@_dry_run
def nbr_matrix(
    paths_cfg, plot_cfg, figpath, force_overwrite, dry_run,
):
    """draws square matrix of sequence logo's from neighbour analysis"""
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

    fig.write(figpath)
    LOGGER.output_file(figpath)
    click.secho("Wrote %s" % figpath, fg="green")


@main.command()
@_fig_cfg
@_figpath
def grid(fig_config, figpath):
    """draws an arbitrary shaped grid of mutation motifs based on a config file"""
    # we read in the config file and determine number of rows and columns
    # paths, headings, etc ..
    # then create the figure and axes and call the mutation_motif drawing code
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
    fig.write(path=figpath)
    click.secho("Wrote Cogent3 %s" % figpath, fg="green")


@main.command()
@_json_path
@_group_label
@_plot_cfg
@_figpath
@_force_overwrite
@_dry_run
def spectra_grid(
    json_path, group_label, plot_cfg, figpath, force_overwrite, dry_run,
):
    """draws logo from mutation spectra analysis"""
    # the following is for logging
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
        json_path, plot_cfg=plot_cfg, group_label=group_label
    )
    fig.write(figpath)
    LOGGER.output_file(figpath)
    click.secho("Wrote %s" % figpath, fg="green")


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
def nbr(
    json_paths, plot_cfg, group_label, force_overwrite, dry_run,
):
    """makes motifs for independent or higher order interactions"""
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
        fig.write(paths[order].outpath)
        LOGGER.output_file(paths[order].outpath)
        click.secho(f"Wrote {paths[order].outpath}", fg="green")

    summary = "summary"
    if summary in paths:
        LOGGER.input_file(paths[summary].inpath)
        fig = get_summary_drawable(paths[summary].inpath, plot_cfg)
        fig.write(paths[summary].outpath)
        LOGGER.output_file(paths[summary].outpath)
        click.secho(f"Wrote {paths[summary].outpath}", fg="green")

    click.secho(f"Done!", fg="green")


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
def mi(
    json_path, plot_cfg, group_label, figpath, use_freq, force_overwrite, dry_run,
):
    """draws conventional sequence logo, using MI, from first order effects"""
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
        data, plot_cfg, group_label=group_label, get_heights=get_mi_plot_data
    )
    for ann in fig.layout.annotations:
        if "RE" in ann.text:
            ann.text = ann.text.replace("RE", "MI")
            break
    fig.write(figpath)
    LOGGER.output_file(figpath)
    click.secho("Wrote %s" % figpath, fg="green")


@main.command()
@click.argument("outpath")
def export_cfg(outpath):
    """exports the sample config files to the nominated path"""
    import shutil

    if os.path.exists(outpath):
        click.secho(
            "outpath already exists, delete it or choose different dest", fg="red"
        )
        sys.exit(1)

    path = resource_filename("mutation_motif", f"cfgs")
    shutil.copytree(path, outpath)
    click.secho("Contents written to %s" % outpath, fg="green")
