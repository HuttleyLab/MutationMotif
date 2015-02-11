#!/usr/bin/env python
import os, glob, re
from itertools import permutations
from optparse import make_option

import numpy

from cogent import LoadSeqs, DNA, LoadTable
from cogent.core.alignment import DenseAlignment
from cogent.util.option_parsing import parse_command_line_parameters

from mutation_motif import profile, util, height, logo, entropy, distribution,\
            motif_count

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = \
"log-linear analysis of neighbouring bases."\
+ " Writes counts tables, estimated statistics and figures to the specified"\
+ "directory outpath."

script_info['required_options'] = [
     make_option('-s','--alignfile', help='fasta aligned file centred on mutated position.'),
     make_option('-o','--outpath', help='Directory path to write data.'),
    make_option('--direction', default=None,
     choices=['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT', 'GtoA', 'GtoC',
             'GtoT', 'TtoA', 'TtoC', 'TtoG'], help='Mutation direction.'),
    ]

script_info['optional_options'] = [
    make_option('--flank_size', type=int, default=5,
        help='Flank size, must be 1 <= flank_size <= 5.'),
    make_option('--coding', action='store_true', default=False,
        help='If True, samples pseudo-SNPs at 3-bp intervals from real SNP.'),
    make_option('--format', default='pdf', choices=['pdf', 'png'],
        help='Plot format.'),
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
    
    assert opts.flank_size <= 5
    
    outpath = util.abspath(opts.outpath)
    if not opts.dry_run:
        util.create_path(outpath)
    
    direction = tuple(opts.direction.split('to'))
    chosen_base = direction[0]
    
    fn = util.abspath(opts.alignfile)
    orig_seqs = util.load_from_fasta(fn)
    seqs = orig_seqs.ArraySeqs
    seqs = util.just_nucs(seqs)
    orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base, step=1,
                                    flank_size=opts.flank_size)
    
    # convert profiles to a motif count table
    orig_counts = motif_count.profile_to_seq_counts(orig, flank_size=opts.flank_size)
    ctl_counts = motif_count.profile_to_seq_counts(ctl, flank_size=opts.flank_size)
    counts_table = motif_count.get_count_table(orig_counts, ctl_counts, opts.flank_size*2)
    
    if not opts.dry_run:
        outpath = os.path.join(outpath, "counts_table.txt")
        counts_table.writeToFile(outpath, sep='\t')
    
    if opts.dry_run or opts.verbose:
        print
        print counts_table
        print
    
    
