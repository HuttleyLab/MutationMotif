"""combines counts from each mutation direction into a single table"""
from __future__ import division

import os, sys, time, re, glob
from collections import Counter

import click

from cogent import LoadTable
from scitrack import CachingLogger

from mutation_motif.util import open_, create_path, abspath, get_subtables
from mutation_motif.complement import make_strand_symmetric_table

LOGGER = CachingLogger(create_dir=True)
_directions = ['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT', 'GtoA', 'GtoC',
                          'GtoT', 'TtoA', 'TtoC', 'TtoG']
direction = re.compile('(%s)' % '|'.join(_directions))

def check_found_filenames(filenames):
    '''check the number of filenames and that they include the direction'''
    found = Counter()
    for fn in filenames:
        d = direction.findall(fn)
        found.update(d)
    total = sum(found.values())
    if total != 12 or set(found) != set(_directions):
        print "ERROR: counts_pattern did not identify 12 files -- %s" % filenames
        print "Note that each file must contain a single direction pattern, e.g. CtoT, AtoG"
        exit(-1)

@click.command()
@click.option('-c','--counts_pattern', help='glob pattern uniquely identifying all 12 mutation counts files.')
@click.option('-o','--output_path', help='Path to write combined_counts data.')
@click.option('-s','--strand_symmetric', is_flag=True,
         help='produces table suitable for strand symmetry test.')
@click.option('-p','--split_dir', default=None,
         help='path to write individual direction strand symmetric tables.')
@click.option('-D','--dry_run', is_flag=True,
        help='Do a dry run of the analysis without writing output.')
@click.option('-F','--force_overwrite', is_flag=True,
        help='Overwrite output and run.log files.')
def main(counts_pattern, output_path, strand_symmetric, split_dir, dry_run, force_overwrite):
    """export tab delimited combined counts table by appending the 12 mutation
    direction tables, adding a new column ``direction``."""
    
    output_path = abspath(output_path)
    if strand_symmetric and split_dir:
        split_dir = abspath(split_dir)
    else:
        split_dir = None
    
    # check we the glob pattern produces the correct number of files
    counts_files = glob.glob(counts_pattern)
    check_found_filenames(counts_files)
    
    counts_filename = os.path.join(output_path, 'combined_counts.txt')
    runlog_path = os.path.join(output_path, "combined_counts.log")
    
    if not dry_run:
        if not force_overwrite and (os.path.exists(counts_filename) or os.path.exists(runlog_path)):
            msg = "Either %s or %s already exist. Force overwrite of existing files with -F."
            raise ValueError(msg % (counts_filename, runlog_path))
        
        create_path(output_path)
        if split_dir:
            create_path(split_dir)
        
        LOGGER.log_file_path = runlog_path
        LOGGER.log_message(str(vars(opts)), label='vars')
        for fn in counts_files:
            LOGGER.input_file(fn, label="count_file")
    
    start_time = time.time()
    
    # run the program
    all_counts = []
    header = None
    num_rows = 0
    basenames = []
    for fn in counts_files:
        basenames.append(os.path.basename(fn))
        mutation = direction.findall(fn)[0]
        table = LoadTable(fn, sep='\t')
        if header is None:
            header = list(table.Header)
            header.append('direction')
            num_rows = table.Shape[0]
        
        data = table.getRawData()
        new = []
        for row in data:
            row.append(mutation)
            new.append(row)
        all_counts += new
    
    table = LoadTable(header=header, rows=all_counts)
    
    if strand_symmetric:
        table = make_strand_symmetric_table(table)
    
    if split_dir:
        group_subtables = get_subtables(table, group_label='direction')
    
    if not dry_run:
        table.writeToFile(counts_filename, sep='\t')
        LOGGER.output_file(counts_filename)
        
        if split_dir:
            for group, subtable in group_subtables:
                # we first assume that group is part of the filenames!
                fn = [bn for bn in basenames if group in bn]
                if len(fn) == 1:
                    fn = fn[0]
                else:
                    fn = "%s.txt" % group
                
                counts_filename = os.path.join(split_dir, fn)
                subtable.writeToFile(counts_filename, sep='\t')
                LOGGER.output_file(counts_filename)
    
    # determine runtime
    duration = time.time() - start_time
    if not dry_run:
        LOGGER.log_message("%.2f" % (duration/60.), label="run duration (minutes)")
    
    print "Done!"
