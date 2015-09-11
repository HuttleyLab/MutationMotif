#!/usr/bin/env python
import os, time, json
from itertools import permutations, combinations
from optparse import make_option

import numpy
from pandas import read_json
from matplotlib import pyplot

from cogent import LoadSeqs, DNA, LoadTable
from cogent.core.alignment import DenseAlignment
from cogent.util.option_parsing import parse_command_line_parameters
from cogent.maths.stats import chisqprob

from scitrack import CachingLogger

from mutation_motif import profile, util, logo, motif_count, log_lin

LOGGER = CachingLogger(create_dir=True)

def dump_json(data, outfile_path):
    with open(outfile_path, mode="w") as outfile:
        json.dump(data, outfile)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "\n".join([
"log-linear analysis of point mutation spectra",
"",
"The spectra is simply the total number of mutations for each point mutation.",
"These counts are compared between groups. Writes estimated",
"statistics, figures and a run log to the specified directory outpath.",
"",
"See documentation for count table format requirements."])

script_info['required_options'] = [
     make_option('-1','--countsfile', help='tab delimited file of counts.'),
     make_option('-o','--outpath', help='Directory path to write data.'),
    ]

script_info['optional_options'] = [
    make_option('-2','--countsfile2',
        help='second group motif counts file.'),
    make_option('-s','--strand_symmetry', action='store_true', default=False,
        help='Second group is strand symmetry.'),
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
    
    table = LoadTable(opts.countsfile, sep='\t')
    if not opts.dry_run:
        log_file_path = os.path.join(util.abspath(opts.outpath),
                                     'spectra_analysis.log')
        LOGGER.log_file_path = log_file_path
    
    LOGGER.input_file(opts.countsfile)
    # if there's a strand symmetry argument then we don't need a second file
    if opts.strand_symmetry:
        group_label = 'strand'
        counts_table = util.spectra_table(table, group_label)
    
    if not opts.strand_symmetry:
        group_label = 'group'
        
        # be sure there's two files
        counts_table2 = LoadTable(opts.countsfile2, sep='\t')
        LOGGER.input_file(opts.countsfile2)
        counts_table2 = counts_table2.withNewColumn('group',
                                lambda x: '2', columns=counts_table2.Header[0])
        counts_table1 = table.withNewColumn('group',
                                lambda x: '1', columns=table.Header[0])
        
        counts_table1 = util.spectra_table(counts_table1, group_label)
        counts_table2 = util.spectra_table(counts_table2, group_label)
        
        # now combine
        header = ['group'] + counts_table2.Header[:-1]
        raw1 = counts_table1.getRawData(header)
        raw2 = counts_table2.getRawData(header)
        counts_table = LoadTable(header=header, rows=raw1+raw2)
        
        if opts.verbose:
            print counts_table
        
    total_re, dev, df, collated, formula = log_lin.spectra_difference(counts_table, group_label)
    r = [list(x) for x in collated.to_records(index=False)]
    if not opts.strand_symmetry:
        grp_labels = {'1': opts.countsfile,
                      '2': opts.countsfile2}
        grp_index = list(collated.columns).index('group')
        for row in r:
            row[grp_index] = grp_labels[row[grp_index]]
    
    
    p = chisqprob(dev, df)
    if p < 1e-6:
        prob = "%.2e" % p
    else:
        prob = "%.6f" % p
    
    significance = ["RE=%.6f" % total_re, "Dev=%.2f" % dev,  "df=%d" % df,
                    "p=%s" % p]
    
    print "  :  ".join(significance)
    
    table = LoadTable(header=collated.columns, rows=r, digits=5).sorted(columns='ret')
    saveable = dict(rel_entropy=total_re, deviance=dev, df=df, prob=p,
                    formula=formula, stats=collated.to_json())
    json_path = None
    
    outpath = util.abspath(opts.outpath)
    if not opts.dry_run:
        util.create_path(outpath)
        json_path = os.path.join(outpath, 'spectra_analysis.json')
        dump_json(saveable, json_path)
        LOGGER.output_file(json_path)
        table_path = os.path.join(outpath, 'spectra_summary.txt')
        table.writeToFile(table_path, sep='\t')
        LOGGER.output_file(table_path)
        LOGGER.write(str(significance), label="significance")


if __name__ == "__main__":
    main()
