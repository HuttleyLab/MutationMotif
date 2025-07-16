import contextlib
import io
import json
import os
import re
import time
from configparser import ConfigParser, NoOptionError, NoSectionError

# to be used as a decorator for click commands
from importlib import resources

import numpy
from cogent3 import DNA, load_table, make_table, open_
from cogent3.util.union_dict import UnionDict
from pandas import read_json


def load_table_from_delimited_file(path, sep="\t"):
    """returns a Table object after a quicker loading"""
    with open_(path, "rt") as infile:
        data = infile.read()

    data = data.splitlines()
    header = data.pop(0).strip().split(sep)
    count_index = header.index("count")
    records = []
    for line in data:
        line = line.strip().split(sep)
        line[count_index] = int(line[count_index])
        records.append(line)
    table = make_table(header=header, rows=records)
    if "direction" in table.columns:
        table = make_consistent_direction_style(table)
    return table


def make_consistent_direction_style(table):
    """returns a table with the XtoY style"""
    current = table.distinct_values("direction")
    pattern = re.compile("[ACGT]to[ACGT]")
    direction_map = {}
    for d in current:
        if "to" in d:
            expect = 4
            val = d
        elif ">" in d:
            expect = 3
            val = f"{d[0]}to{d[2]}"
        else:
            msg = f"unknown direction '{d}'"
            raise ValueError(msg)

        assert len(d) == expect, f"unknown direction '{d}'"
        match = pattern.search(val)
        if match is None:
            msg = f"unknown direction '{d}'"
            raise ValueError(msg)

        direction_map[d] = val

    new_col = [direction_map[v] for v in table.columns["direction"]]
    table.columns["direction"] = numpy.array(new_col, dtype="U")
    return table


def spectra_table(table, group_label):
    """returns a table with columns without position information"""
    assert "direction" in table.header
    if "mut" in table.header:
        # remove redundant category (counts of M == R)
        table = table.filtered("mut=='M'")

    columns = ["count", "direction", group_label]
    table = table.get_columns(columns)
    # so we have a table with counts per direction
    results = []
    group_categories = table.distinct_values(group_label)
    filter_template = "direction=='%(direction)s' and %(label)s=='%(category)s'"
    for direction in table.distinct_values("direction"):
        start = direction[0]
        for group_category in group_categories:
            condition = {
                "direction": direction,
                "label": group_label,
                "category": group_category,
            }
            sub_table = table.filtered(filter_template % condition)
            total = sub_table.summed("count")
            results.append([total, start, direction, group_category])
    result = make_table(
        header=["count", "start", "direction", group_label],
        rows=results,
    )
    return make_consistent_direction_style(result)


def get_subtables(table, group_label="direction"):
    """returns [(group, subtable),...] for distinct values of group_label"""
    groups = table.distinct_values(group_label)
    tables = []
    for group in groups:
        subtable = table.filtered(lambda x: x == group, columns=group_label)
        tables.append((group, subtable))
    return tables


def dump_loglin_stats(data, outfile_path) -> None:
    """save data in json format to outfile_path"""
    # convert all pandas data frames to json
    saveable = {}
    for position_set in data:
        curr = {}
        new_key = str(position_set)
        for key, value in list(data[position_set].items()):
            if key == "stats":
                value = data[position_set][key].to_json()
            curr[key] = value

        saveable[new_key] = curr

    with open(outfile_path, mode="w") as outfile:
        json.dump(saveable, outfile)


def load_loglin_stats(infile_path):
    """read in data in json format"""
    # convert all 'stats' to pandas data frames
    with open(infile_path) as infile:
        data = json.load(infile)

    new_data = {}
    for position_set in data:
        try:
            new_key = eval(position_set)
        except NameError:
            new_key = position_set

        new_data[new_key] = {}
        for key, value in list(data[position_set].items()):
            if key == "stats":
                with io.StringIO(value) as stream:
                    value = read_json(stream)
            new_data[new_key][key] = value
    return new_data


def is_valid(data):
    """returns True if all elements satisfy 0 <= e < 4"""
    return (data >= 0).all() and (data < 4).all()


def array_to_str(data):
    """convert numpy array back to DNA sequence"""
    data = data.astype(numpy.uint8)
    return ["".join(DNA.alphabet.from_indices(v)) for v in data]


def seqs_to_array(d_aln):
    """get input array alignment data, transfer them into a numpy array matrix,
    whith accordant numbers for DNA bases
    also
    filter sequences and save just_nuc sequences.
    """
    return just_nucs(d_aln.array_seqs)


def just_nucs(seqs):
    """eliminate sequences containing gaps/Ns,
    along the 1st axis, match each base in seqs to <= 3 element-wise,
    give the indices of those just_nucs seq idx.
    """
    (indices,) = (seqs <= 3).all(axis=1).nonzero()
    return seqs.take(indices, axis=0)


def abspath(path):
    """returns an expanded, absolute path"""
    return os.path.abspath(os.path.expanduser(path))


def makedirs(path) -> None:
    """creates dir path"""
    with contextlib.suppress(OSError):
        os.makedirs(path, exist_ok=True)


def get_selected_indices(stats, group_label=None, group_ref=None):
    """returns indices for selecting dataframe records for display"""
    if group_label and group_ref is None:  # TODO this logic needs improving
        val = {"strand": "+"}.get(group_label, "1")
        return numpy.logical_and(stats["mut"] == "M", stats[group_label] == val)
    if group_label and group_ref:
        return numpy.logical_and(
            stats["mut"] == "M",
            stats[group_label] == group_ref,
        )
    return stats["mut"] == "M"


_pos_num = re.compile(r"\d+$")


def get_position_number(pos):
    """returns the position index from a string formatted as 'pos##'"""
    return int(_pos_num.search(pos).group())


def get_order_max_re_from_summary(table):
    """returns the maximum RE for each effect order"""
    if isinstance(table, str):
        table = load_table(table, sep="\t")
    table = table.with_new_column(
        "order",
        lambda x: x.count(":") + 1,
        columns="Position",
    )
    orders = table.distinct_values("order")
    table = table.get_columns(["order", "RE"])
    rows = []
    for order in orders:
        subtable = table.filtered(lambda x: x == order, columns="order")
        rows.append([order, max(subtable.to_list("RE"))])
    return rows


def get_config_parser(path, default):
    if not path or not os.path.exists(path):
        path = resources.files("mutation_motif") / f"cfgs/{default}"

    parser = ConfigParser()
    parser.optionxform = str  # stops automatic conversion to lower case
    if isinstance(path, io.TextIOBase):
        parser.read_file(path)
    else:
        parser.read(str(path))

    return parser


def get_fig_properties(parser, section="fig setup"):
    get_val = parser.get
    cfg = UnionDict()
    figsize = [
        100 * float(v.strip()) / 2.5 for v in get_val(section, "figsize").split(",")
    ]
    cfg.width, cfg.height = figsize

    with contextlib.suppress(NoOptionError):
        margin = {
            k[0]: int(get_val(section, k)) for k in ("top", "bottom", "right", "left")
        }
        cfg.margin = margin

    # font sizes, title, label text padding
    for option in parser.options(section):
        valid = any(attr in option for attr in ("pad", "font", "angle"))
        if not valid:
            continue

        cfg[option] = int(get_val(section, option))

    # ylim
    with contextlib.suppress(NoOptionError):
        ylim = float(get_val(section, "ylim"))
        cfg.ylim = ylim

    with contextlib.suppress(NoOptionError):
        ylabel = get_val(section, "ylabel")
        cfg.ylabel = ylabel

    with contextlib.suppress(NoOptionError):
        space = float(get_val(section, "space"))
        cfg.space = space

    with contextlib.suppress(NoOptionError):
        xlabel = get_val(section, "xlabel")
        cfg.xlabel = xlabel

    for axis in ("rows", "cols"):
        key = f"num_{axis}"
        with contextlib.suppress(NoOptionError):
            val = int(parser.get(section, key))
            cfg[key] = val

    return cfg


def get_spectra_config(path):
    parser = get_config_parser(path, default="spectra.cfg")
    col_title = parser.get("fig setup", "col_title")
    row_title = parser.get("fig setup", "row_title")
    cfg = get_fig_properties(parser, section="fig setup")
    cfg.col_title = col_title
    cfg.row_title = row_title
    return cfg


def get_nbr_config(path, section):
    parser = get_config_parser(path, default="nbr.cfg")
    return get_fig_properties(parser, section=section)


def get_nbr_matrix_config(path):
    parser = get_config_parser(path, default="nbr_matrix.cfg")
    cfg = get_fig_properties(parser, section="fig setup")
    col_title = parser.get("fig setup", "col_title")
    row_title = parser.get("fig setup", "row_title")
    cfg.col_title = col_title
    cfg.row_title = row_title

    return cfg


def get_grid_config(path):
    dirname = None if path is None else os.path.dirname(path)
    parser = get_config_parser(path, default="grid.cfg")
    cfg = get_fig_properties(parser, section="fig setup")
    subplots = UnionDict()
    if path:
        # load user defined values
        with contextlib.suppress(NoOptionError):
            col_titles = [
                l.strip() for l in parser.get("fig setup", "col_titles").split(",")
            ]
            cfg.col_titles = col_titles
            num_titles = len(col_titles)
            if num_titles != cfg.num_cols:
                msg = f"number of col_titles {num_titles} != num_cols {cfg.num_cols}"
                raise ValueError(msg)

        with contextlib.suppress(NoOptionError):
            row_titles = [
                l.strip() for l in parser.get("fig setup", "row_titles").split(",")
            ]
            cfg.row_titles = row_titles
            num_titles = len(row_titles)
            if num_titles != cfg.num_rows:
                msg = f"number of row_titles {num_titles} != num_rows {cfg.num_rows}"
                raise ValueError(msg)

        # load possible sections
        for section in parser.sections():
            try:
                coord = tuple(i - 1 for i in map(int, section.split(",")))
                if coord[0] >= (cfg.num_cols):
                    msg = f"subplot section exceeds specified num_cols={cfg.num_cols}"
                    raise ValueError(msg)
                if coord[1] >= (cfg.num_rows):
                    msg = f"subplot section exceeds specified num_rows={cfg.num_rows}"
                    raise ValueError(msg)
            except (TypeError, ValueError):
                continue

            try:
                path = parser.get(section, "path")
            except NoOptionError:
                continue

            if path is None:
                continue

            if dirname and dirname not in path:
                path = os.path.join(dirname, path)
            subplots[coord] = path

    if not subplots and path is not None:
        msg = "no sections for subplots"
        raise NoSectionError(msg)

    cfg.paths = subplots

    return cfg


def get_summary_config(path):
    parser = get_config_parser(path, default="nbr.cfg")
    return get_fig_properties(parser, section="summary")


def get_nbr_path_config(path):
    dirname = None if path is None else os.path.dirname(path)
    parser = get_config_parser(path, default="nbr_paths.cfg")
    cfg = UnionDict()
    for section in parser.sections():
        paths = {k: parser.get(section, k) for k in ("inpath", "outpath")}
        if not paths.get("inpath"):
            continue

        if dirname and dirname not in paths["inpath"]:
            inpath = os.path.join(dirname, paths["inpath"])
            paths["inpath"] = inpath
            if section != "summary":
                assert inpath.endswith(".json"), f"{inpath} missing json suffix"

        outpath = paths.get("outpath")
        if not outpath:
            outpath = paths["inpath"].replace(".json", ".pdf")
        elif dirname and dirname not in outpath:
            outpath = os.path.join(dirname, outpath)
        paths["outpath"] = outpath

        cfg[section] = UnionDict(paths)
    return cfg


def est_ylim(char_heights):
    """returns a ylim for character height axis plotting"""
    try:
        t = numpy.fabs(char_heights).sum(axis=1).max()
    except ValueError:
        t = numpy.fabs(char_heights).max()

    for i in range(1, 10):
        ylim = numpy.around(t, i)
        if ylim != 0:
            break

    if ylim < t:
        ylim *= 1.667
        ylim = numpy.around(ylim, i)

    return max(ylim, 1e-6)


class pdf_writer:
    """class that handles super annooying mathjax warning box in plotly pdf's"""

    def __init__(self) -> None:
        self._done_once = False

    def __call__(self, fig, path):
        # the sleep, plus successive write, is ESSENTIAL to avoid the super annoying
        # "[MathJax]/extensions/MathMenu.js" text box error
        # but we only need to do this once
        if not self._done_once:
            fig.write(path)
            time.sleep(2)
            self._done_once = True
        fig.write(path)
