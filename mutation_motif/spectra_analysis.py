#!/usr/bin/env python
import os
import json

from cogent3 import LoadTable
from cogent3.maths.stats import chisqprob

from scitrack import CachingLogger

from mutation_motif import util, log_lin

LOGGER = CachingLogger(create_dir=True)


def dump_json(data, outfile_path):
    with open(outfile_path, mode="w") as outfile:
        json.dump(data, outfile)


def main(countsfile, outpath, countsfile2, strand_symmetry, force_overwrite,
         dry_run, verbose):
    args = locals()

    table = LoadTable(countsfile, sep='\t')
    if not dry_run:
        log_file_path = os.path.join(util.abspath(outpath),
                                     'spectra_analysis.log')
        LOGGER.log_file_path = log_file_path
        LOGGER.log_message(str(args), label='vars')

    LOGGER.input_file(countsfile)
    # if there's a strand symmetry argument then we don't need a second file
    if strand_symmetry:
        group_label = 'strand'
        counts_table = util.spectra_table(table, group_label)

    if not strand_symmetry:
        group_label = 'group'

        # be sure there's two files
        counts_table2 = LoadTable(countsfile2, sep='\t')
        LOGGER.input_file(countsfile2)
        counts_table2 = counts_table2.with_new_column('group',
                                                      lambda x: '2', columns=counts_table2.header[0])
        counts_table1 = table.with_new_column('group',
                                              lambda x: '1', columns=table.header[0])

        counts_table1 = util.spectra_table(counts_table1, group_label)
        counts_table2 = util.spectra_table(counts_table2, group_label)

        # now combine
        header = ['group'] + counts_table2.header[:-1]
        raw1 = counts_table1.tolist(header)
        raw2 = counts_table2.tolist(header)
        counts_table = LoadTable(header=header, rows=raw1 + raw2)

        if verbose:
            print(counts_table)

    # spectra table has [count, start, end, group] order
    # we reduce comparisons to a start base
    results = []
    saveable = {}
    for start_base in counts_table.distinct_values('start'):
        subtable = counts_table.filtered('start == "%s"' % start_base)
        columns = [c for c in counts_table.header if c != 'start']
        subtable = subtable.get_columns(columns)
        total_re, dev, df, collated, formula = log_lin.spectra_difference(
            subtable, group_label)
        r = [list(x) for x in collated.to_records(index=False)]

        if not strand_symmetry:
            grp_labels = {'1': countsfile,
                          '2': countsfile2}
            grp_index = list(collated.columns).index('group')
            for row in r:
                row[grp_index] = grp_labels[row[grp_index]]

        p = chisqprob(dev, df)
        if p < 1e-6:
            prob = "%.2e" % p
        else:
            prob = "%.6f" % p

        for row in r:
            row.insert(0, start_base)
            row.append(prob)

        results += r

        significance = ["RE=%.6f" % total_re, "Dev=%.2f" % dev, "df=%d" % df,
                        "p=%s" % p]

        stats = "  :  ".join(significance)
        print("Start base=%s  %s" % (start_base, stats))
        saveable[start_base] = dict(rel_entropy=total_re, deviance=dev,
                                    df=df, prob=p,
                                    formula=formula, stats=collated.to_json())

    table = LoadTable(header=['start_base'] + list(collated.columns) +
                             ['prob'],
                      rows=results, digits=5).sorted(columns='ret')
    json_path = None

    outpath = util.abspath(outpath)
    if not dry_run:
        util.create_path(outpath)
        json_path = os.path.join(outpath, 'spectra_analysis.json')
        dump_json(saveable, json_path)
        LOGGER.output_file(json_path)
        table_path = os.path.join(outpath, 'spectra_summary.txt')
        table.write(table_path, sep='\t')
        LOGGER.output_file(table_path)
        LOGGER.log_message(str(significance), label="significance")
