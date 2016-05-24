"""combines counts from each mutation direction into a single table"""
from __future__ import division

import os, sys, time, re, glob
from collections import Counter
from optparse import make_option

from cogent import LoadTable
from cogent.util.option_parsing import parse_command_line_parameters
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
    

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "export tab delimited combined counts "\
"table by appending the 12 mutation direction tables, adding a new column "\
"``direction``."

script_info['required_options'] = [
     make_option('-c','--counts_pattern', help='glob pattern uniquely identifying all 12 mutation counts files.'),
     make_option('-o','--output_path', help='Path to write combined_counts data.'),
    ]

script_info['optional_options'] = [
    make_option('-s','--strand_symmetric', action='store_true', default=False,
         help='produces table suitable for strand symmetry test.'),
    make_option('-p','--split_dir', default=None,
         help='path to write individual direction strand symmetric tables.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    make_option('-F','--force_overwrite', action='store_true', default=False,
        help='Overwrite output and run.log files.'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    opts.output_path = abspath(opts.output_path)
    if opts.strand_symmetric and opts.split_dir:
        split_dir = abspath(opts.split_dir)
    else:
        split_dir = None
    
    # check we the glob pattern produces the correct number of files
    counts_files = glob.glob(opts.counts_pattern)
    check_found_filenames(counts_files)
    
    counts_filename = os.path.join(opts.output_path, 'combined_counts.txt')
    runlog_path = os.path.join(opts.output_path, "combined_counts.log")
    
    if not opts.dry_run:
        if not opts.force_overwrite and (os.path.exists(counts_filename) or os.path.exists(runlog_path)):
            msg = "Either %s or %s already exist. Force overwrite of existing files with -F."
            raise ValueError(msg % (counts_filename, runlog_path))
        
        create_path(opts.output_path)
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
    
    if opts.strand_symmetric:
        table = make_strand_symmetric_table(table)
    
    if split_dir:
        group_subtables = get_subtables(table, group_label='direction')
    
    if not opts.dry_run:
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
    if not opts.dry_run:
        LOGGER.log_message("%.2f" % (duration/60.), label="run duration (minutes)")
    
    print "Done!"
