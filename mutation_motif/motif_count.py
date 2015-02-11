from collections import Counter
from itertools import product

from cogent import DNA, LoadTable

from mutation_motif.util import array_to_str

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
