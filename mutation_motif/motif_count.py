from collections import Counter
from itertools import product

from cogent import DNA, LoadTable

from mutation_motif.util import array_to_str, load_from_fasta, just_nucs
from mutation_motif.profile import get_profiles

def counts_from_seqs(fn, chosen_base, flank_size, seed):
    """returns a counts table of motifs for mutated, control seqs"""
    orig_seqs = load_from_fasta(fn)
    seqs = orig_seqs.ArraySeqs
    seqs = just_nucs(seqs)
    orig, ctl = get_profiles(seqs, chosen_base=chosen_base, step=1,
                                    flank_size=flank_size, seed=seed)

    # convert profiles to a motif count table
    orig_counts = profile_to_seq_counts(orig, flank_size=flank_size)
    ctl_counts = profile_to_seq_counts(ctl, flank_size=flank_size)
    counts_table = get_count_table(orig_counts, ctl_counts,
                                                flank_size*2)
    counts_table = counts_table.sorted(columns='mut')
    return counts_table

def profile_to_seq_counts(data, flank_size):
    """converts data to seqs and returns sequence counts"""
    mp = data.shape[1] // 2
    indices = [i for i in range(mp-flank_size, mp+flank_size+1) if i != mp]
    
    data = data.take(indices, axis=1)
    seqs = array_to_str(data)
    seqs = Counter(seqs)
    return seqs

def get_count_table(observed, control, k):
    """return table of motif counts
    
    Each motif position is a separate column. All possible DNA motifs of length
    k are included.
    
    Arguments:
        - observed: the observed counts as {seq: count}
        - control: the control counts as {seq: count}
        - k: size of the motif"""
    rows = []
    for state in product('ACGT', repeat=k):
        state = ''.join(state)
        control_counts = control[state]
        observed_counts = observed[state]
        rows.append([control_counts] + list(state) + ['U'])
        rows.append([observed_counts] + list(state) + ['M'])
    
    header=['count'] + ["pos%d" % i for i in range(k)] + ['mut']
    table = LoadTable(header=header, rows=rows)
    return table

def reduced_multiple_positions(table, *positions):
    counts = Counter()
    for row in table:
        count = row['count']
        motif = tuple(row[pos] for pos in positions)
        counts[motif] += count
    return counts

def reduced_one_position(table, pos):
    '''returns base counts for one position'''
    base_counts = Counter()
    for row in table:
        count = row['count']
        base = row[pos]
        base_counts[base] += count
    return base_counts

def get_combined_counts(table, positions):
    mutated = table.filtered("mut=='M'")
    unmutated = table.filtered("mut=='U'")
    if type(positions) == str:
        mut_counts = reduced_one_position(mutated, positions)
        unmut_counts = reduced_one_position(unmutated, positions)
        positions = [positions]
        states = 'ACGT'
        header = ['mut', 'base', 'count']
    else:
        mut_counts = reduced_multiple_positions(mutated, *positions)
        unmut_counts = reduced_multiple_positions(unmutated, *positions)
        states = product('ACGT', repeat=len(positions))
        header = ['mut'] + ['base%d' % (i+1) for i in range(len(positions))] + ['count']
        
    
    combined = []
    n = 0
    for state in states:
        n += 1
        combined.append(['U'] + list(state) + [unmut_counts[state]])
        combined.append(['M'] + list(state) + [mut_counts[state]])
    
    counts_table = LoadTable(header=header, rows=combined)
    counts_table = counts_table.sorted(columns=header[:-1])
    return counts_table

