import bz2
import gzip
import io
import json
import os
import re
from configparser import (ConfigParser, NoOptionError, NoSectionError,
                          ParsingError, RawConfigParser)

import click
import numpy
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
from numpy import around
from numpy.core._multiarray_umath import fabs
from pandas import read_json
# to be used as a decorator for click commands
from pkg_resources import resource_filename

from cogent3 import DNA, load_table, make_table
from cogent3.core.alignment import ArrayAlignment
from cogent3.parse.fasta import MinimalFastaParser

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

# to be used as a decorator for click commands
no_type3_font = click.option('--no_type3', is_flag=True,
                             help='Exclude Type 3 fonts from pdf, necessary'
                             ' for ScholarOne figures')


def exclude_type3_fonts():
    """stops matplotlib from using Type 3 fonts"""
    # this is because ScholarOne does not allow embedding of Type 3 fonts in
    # pdf's, a royal PITA. Thanks ScholarOne!!

    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42


def get_plot_configs(cfg_path=None):
    """returns a config object with plotting settings"""
    defaults = dict(xlabel_fontsize=14, ylabel_fontsize=14,
                    xtick_fontsize=12, ytick_fontsize=12,
                    xlabel_pad=0.01, ylabel_pad=0.01)

    figwidths = {'1-way plot': 2.25, '2-way plot': 9, '3-way plot': 9,
                 '4-way plot': 9, 'summary plot': 9, 'grid': 8}
    figsizes = {'1-way plot': (9, 3), '2-way plot': (9, 9),
                '3-way plot': (9, 9), '4-way plot': (9, 9),
                'summary plot': (9, 9),
                'grid': (8, 8)}

    config = RawConfigParser()
    for section in ['1-way plot', '2-way plot', '3-way plot', '4-way plot',
                    'summary plot', 'grid']:
        config.add_section(section)
        config.set(section, 'figwidth', figwidths[section])
        config.set(section, 'figsize', figsizes[section])
        for arg, default in list(defaults.items()):
            config.set(section, arg, default)
    if cfg_path:
        # load the user provided config
        user_config = RawConfigParser(allow_no_value=True)
        try:
            user_config.read(cfg_path)
        except ParsingError as err:
            msg = 'Could not parse %s: %s' % (cfg_path, err)
            raise ParsingError(msg)

        # update the default settings
        for section in config.sections():
            for key, val in config.items(section):
                try:
                    new_val = user_config.get(section, key)
                    config.set(section, key, eval(new_val))
                except (NoSectionError, NoOptionError):
                    pass
    return config


# Code from Joe Kington's answer to
# http://stackoverflow.com/questions/3677368/matplotlib-format-axis-offset-values-to-whole-numbers-or-specific-number
class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""

    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset,
                                 useMathText=useMathText)

    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag


def load_table_from_delimited_file(path, sep='\t'):
    '''returns a Table object after a quicker loading'''
    with open_(path, 'rt') as infile:
        header = infile.readline().strip().split(sep)
        count_index = header.index('count')
        records = []
        for line in infile:
            line = line.strip().split(sep)
            line[count_index] = int(line[count_index])
            records.append(line)
        table = make_table(header=header, rows=records)
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
    filter_template = "direction=='%(direction)s' and " "%(label)s=='%(category)s'"
    for direction in table.distinct_values("direction"):
        start = direction[0]
        for group_category in group_categories:
            condition = dict(
                direction=direction, label=group_label, category=group_category
            )
            sub_table = table.filtered(filter_template % condition)
            total = sub_table.summed("count")
            results.append([total, start, direction, group_category])
    result = make_table(
        header=["count", "start", "direction", group_label], rows=results
    )
    return result


def get_subtables(table, group_label="direction"):
    """returns [(group, subtable),...] for distinct values of group_label"""
    groups = table.distinct_values(group_label)
    tables = []
    for group in groups:
        subtable = table.filtered(lambda x: x == group, columns=group_label)
        tables.append((group, subtable))
    return tables


def dump_loglin_stats(data, outfile_path):
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
                value = read_json(value)
            new_data[new_key][key] = value
    return new_data


def is_valid(data):
    """returns True if all elements satisfy 0 <= e < 4"""
    return (data >= 0).all() and (data < 4).all()


def load_from_fasta(filename):
    infile = open_(filename, mode="rt")
    parser = MinimalFastaParser(infile)
    seqs = [(n, s) for n, s in parser]
    infile.close()
    return ArrayAlignment(data=seqs, moltype=DNA)


def array_to_str(data):
    """convert numpy array back to DNA sequence"""
    return ["".join(DNA.alphabet.from_indices(v)) for v in data]


def seqs_to_array(d_aln):
    """get input array alignment data, transfer them into a numpy array matrix,
    whith accordant numbers for DNA bases
    also
    filter sequences and save just_nuc sequences.
    """
    just_bases = just_nucs(d_aln.array_seqs)
    return just_bases


def just_nucs(seqs):
    """eliminate sequences containing gaps/Ns,
    along the 1st axis, match each base in seqs to <= 3 element-wise,
    give the indices of those just_nucs seq idx.
    """
    (indices,) = (seqs <= 3).all(axis=1).nonzero()
    just_bases = seqs.take(indices, axis=0)
    return just_bases


def open_(filename, mode="r"):
    """handles different compression"""
    op = {"gz": gzip.open, "bz2": bz2.BZ2File}.get(filename.split(".")[-1], open)
    return op(filename, mode)


def abspath(path):
    """returns an expanded, absolute path"""
    return os.path.abspath(os.path.expanduser(path))


def makedirs(path):
    """creates dir path"""
    try:
        os.makedirs(path)
    except OSError as e:
        pass
