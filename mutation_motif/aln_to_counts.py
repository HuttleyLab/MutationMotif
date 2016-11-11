"""export seq files for different mutation types"""


import os
import time
import re

import click

from scitrack import CachingLogger

from mutation_motif.util import makedirs, abspath, just_nucs,\
    load_from_fasta
from mutation_motif import profile, motif_count

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


LOGGER = CachingLogger(create_dir=True)
fn_suffixes = re.compile(r"\.(fa|fasta)\.*(gz|gzip|bz2)*$")


def get_counts_filename(align_path, output_dir):
    """returns counts output path

    Arguments:
        - align_path: path to the alignment file. The basename will be
          modified to use a .txt suffix
        - output_dir: directory where the counts file is to be written
    """

    fn = os.path.basename(align_path)
    fn = fn_suffixes.sub(".txt", fn)
    counts_filename = os.path.join(output_dir, fn)
    return counts_filename


def align_to_counts(align_path, output_path, flank_size, direction, step,
                    seed, randomise, dry_run):
    '''returns counts table from alignment of sequences centred on a SNP'''

    if not dry_run:
        makedirs(output_path)

    print("Deriving counts from sequence file")
    step = int(step)

    direction = tuple(direction.split('to'))
    chosen_base = direction[0]
    orig_seqs = load_from_fasta(os.path.abspath(align_path))
    seqs = orig_seqs.array_seqs
    seqs = just_nucs(seqs)
    if not randomise:
        orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base,
                                         step=step, flank_size=flank_size,
                                         seed=seed)
    else:
        LOGGER.log_message("A randomised selection of starting base "
                           "locations use for observed counts.")
        # we are setting a randomised set of locations as our observed SNPs
        ctl = profile.get_control(seqs, chosen_base=chosen_base, step=step,
                                  flank_size=flank_size, seed=seed)
        orig = profile.get_control(seqs, chosen_base=chosen_base, step=step,
                                   flank_size=flank_size, seed=seed)

    # convert profiles to a motif count table
    orig_counts = motif_count.profile_to_seq_counts(
        orig, flank_size=flank_size)
    ctl_counts = motif_count.profile_to_seq_counts(ctl, flank_size=flank_size)
    counts_table = motif_count.get_count_table(
        orig_counts, ctl_counts, flank_size * 2)
    counts_table = counts_table.sorted(columns='mut')
    return counts_table


@click.command()
@click.option("-a", "--align_path", required=True,
              help="fasta aligned file centred on mutated position.")
@click.option('-o', '--output_path', required=True, help='Path to write data.')
@click.option('-f', '--flank_size', required=True, type=int,
              help='Number of bases per side to include.')
@click.option('--direction', required=True, default=None,
              type=click.Choice(['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG',
                                 'CtoT', 'GtoA', 'GtoC', 'GtoT', 'TtoA',
                                 'TtoC', 'TtoG']),
              help='Mutation direction.')
@click.option('-S', '--seed',
              help='Seed for random number generator (e.g. 17, or 2015-02-13).'
              ' Defaults to system time.')
@click.option('-R', '--randomise', is_flag=True,
              help='Randomises the observed data, observed and reference '
              'counts distributions should match.')
@click.option('--step', default='1', type=click.Choice(['1', '2', '3']),
              help='Specifies a "frame" for selecting the random base.')
@click.option('-D', '--dry_run', is_flag=True,
              help='Do a dry run of the analysis without writing output.')
@click.option('-F', '--force_overwrite', is_flag=True,
              help='Overwrite output and run.log files.')
def main(align_path, output_path, flank_size, direction, seed, randomise,
         step, dry_run, force_overwrite):
    """Export tab delimited counts table from alignment centred on SNP position.

    Output file is written to the same path with just the file suffix changed
    from fasta to txt."""
    args = locals()
    if not seed:
        seed = str(time.time())

    align_path = abspath(align_path)
    output_path = abspath(output_path)

    counts_filename = get_counts_filename(align_path, output_path)
    runlog_path = counts_filename.replace(".txt", ".log")
    LOGGER.log_file_path = runlog_path

    if not dry_run:
        if not force_overwrite and (os.path.exists(counts_filename) or
                                    os.path.exists(runlog_path)):
            msg = "Either %s or %s already exist. Force overwrite of existing"\
                  " files with -F."
            raise ValueError(msg % (counts_filename, runlog_path))

        makedirs(output_path)

        LOGGER.log_message(str(args), label='vars')
        LOGGER.input_file(align_path, label="align_path")
        LOGGER.log_message(str(seed), label="random_number_seed")

    start_time = time.time()

    # run the program

    counts_table = align_to_counts(align_path, output_path, flank_size,
                                   direction, step, seed, randomise, dry_run)

    if not dry_run:
        counts_table.write(counts_filename, sep='\t')
        LOGGER.output_file(counts_filename)

    # determine runtime
    duration = time.time() - start_time
    if not dry_run:
        LOGGER.log_message("%.2f" % (duration / 60.),
                           label="run duration (minutes)")
