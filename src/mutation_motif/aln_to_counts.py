"""export seq files for different mutation types"""

import os
import re

from mutation_motif import motif_count, profile
from mutation_motif.util import just_nucs, load_from_fasta, makedirs

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
    return os.path.join(output_dir, fn)


def align_to_counts(
    align_path,
    output_path,
    flank_size,
    direction,
    step,
    seed,
    randomise,
    dry_run,
    LOGGER,
):
    """returns counts table from alignment of sequences centred on a SNP"""

    if not dry_run:
        makedirs(output_path)

    print("Deriving counts from sequence file")
    step = int(step)

    direction = tuple(direction.split("to"))
    chosen_base = direction[0]
    orig_seqs = load_from_fasta(os.path.abspath(align_path))
    seqs = orig_seqs.array_seqs
    seqs = just_nucs(seqs)
    if not randomise:
        orig, ctl = profile.get_profiles(
            seqs,
            chosen_base=chosen_base,
            step=step,
            flank_size=flank_size,
            seed=seed,
        )
    else:
        LOGGER.log_message(
            "A randomised selection of starting base "
            "locations use for observed counts.",
        )
        # we are setting a randomised set of locations as our observed SNPs
        ctl = profile.get_control(
            seqs,
            chosen_base=chosen_base,
            step=step,
            flank_size=flank_size,
            seed=seed,
        )
        orig = profile.get_control(
            seqs,
            chosen_base=chosen_base,
            step=step,
            flank_size=flank_size,
            seed=seed,
        )

    # convert profiles to a motif count table
    orig_counts = motif_count.profile_to_seq_counts(orig, flank_size=flank_size)
    ctl_counts = motif_count.profile_to_seq_counts(ctl, flank_size=flank_size)
    counts_table = motif_count.get_count_table(orig_counts, ctl_counts, flank_size * 2)
    counts_table = counts_table.sorted(columns="mut")
    return counts_table
