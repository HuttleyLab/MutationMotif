from collections import Counter
from itertools import product
from numpy import logical_and, flatnonzero, zeros, array, matrix, prod, \
    ndarray, asarray
from random import choice, seed as set_seed

from cogent import LoadSeqs, DNA, LoadTable
from cogent.align.weights.util import AlnToProfile
from cogent.core.alignment import DenseAlignment, seqs_from_array

from mutation_motif.util import seqs_to_array as load_seqs_to_array, array_to_str

def get_zero_counts(dim, dtype, pseudo_count=1):
    """returns a ProfileCounts instance zeroed"""
    data = zeros((4, dim), dtype=dtype)
    return ProfileCounts(data, pseudo_count=pseudo_count)

class ProfileCounts(ndarray):
    """counts object"""
    def __new__(cls, data, pseudo_count=0):
        new = asarray(data).view(cls)
        new += pseudo_count
        return new
    
    def add_seq(self, seq):
        """add a new sequence"""
        for i in range(seq.shape[0]):
            self[seq[i], i] += 1
    
    

def MakeCircleRange(circle_size, slice_side):
    """factory function that pre computes all circular slices
    slices are assumed to be 2*slice_side+1 in length
    """
    assert (circle_size - 1) % 2 == 0

    slice_collection = {}
    full_indices = list(range(0, circle_size))
    for centre_index in range(circle_size):
        left = centre_index - slice_side
        right = centre_index + slice_side + 1

        if left < 0:
            indices = full_indices[left:] + full_indices[:right]
        elif right > circle_size:
            indices = full_indices[left:] + full_indices[:right - circle_size]
        else:
            indices = full_indices[left: right]
        slice_collection[centre_index] = indices

    def call(centre_index):
        return slice_collection[centre_index]

    return call

def chosen_base_indices(data, chosen_base, step):
    """returns a list of arrays containing positions with chosen base (starting base)"""
    assert (data.shape[1]-1) % 2 == 0, "seqs not 2n + 1 long"
    if type(chosen_base) == str:
        assert len(chosen_base) == 1, "single base only"

        chosen_base = DNA.Alphabet.toIndices(chosen_base)[0]

    indices = zeros(data.shape, dtype=bool)

    # define in-frame columns
    for column in range(data.shape[1] / 2 % step, data.shape[1], step):
        if column == data.shape[1] / 2:
            continue
        indices[:, column] = True

    indices = logical_and(data == chosen_base, indices)
    indices_by_row = []
    for i in range(data.shape[0]):
        row_indices, = indices[i, :].nonzero()
        indices_by_row.append(row_indices)

    # exclude seqs with minimum number
    return indices_by_row

def filter_seqs_by_chosen_base(data, chosen_base_indices, min_chosen_bases):
    """returns sequences and end base indices for seqs with sufficient end bases"""
    keep = []
    new_indices = []
    for i in range(data.shape[0]):
        if len(chosen_base_indices[i]) >= min_chosen_bases:
            keep.append(i)
            new_indices.append(chosen_base_indices[i])

    data = data.take(keep, axis=0)
    return data, new_indices

def get_random_indices(chosen_base_indices, circle_range):
    sampled = []
    for v in chosen_base_indices:
        r = choice(v)
        sampled.append(circle_range(r))

    return sampled

def MakeControl(data, chosen_base, step, flank_size):
    """factory function for generating control profiles"""
    sample_indices = chosen_base_indices(data, chosen_base, step)
    data, sample_indices = filter_seqs_by_chosen_base(data, sample_indices, 1)
    circle_range = MakeCircleRange(data.shape[1], flank_size)
    def call():
        return get_control(data, chosen_base, step, flank_size,
                sample_indices=sample_indices, circle_range=circle_range)
    return call

def get_control(seq_array, chosen_base, step, flank_size, sample_indices=None, circle_range=None, seed=None):
    assert seed is not None, "Must provide a random number seed"
    set_seed(seed)
    if sample_indices is None:
        sample_indices = chosen_base_indices(seq_array, chosen_base, step)
        seq_array, sample_indices = filter_seqs_by_chosen_base(seq_array, sample_indices, 1)
    
    if circle_range is None:
        circle_range = MakeCircleRange(seq_array.shape[1], flank_size)
    
    sampled_indices = get_random_indices(sample_indices, circle_range)
    rows = []
    for i in range(len(seq_array)):
        row = seq_array[i].take(sampled_indices[i])
        rows.append(row)
    sampled_data = array(rows)
    return sampled_data

def get_profiles(data, chosen_base, step, flank_size, circle_range=None, seed=None):
    """returns matched profile and control profile"""
    ctl = get_control(data, chosen_base, step, flank_size,
                        circle_range=circle_range, seed=seed)
    length = data.shape[1]
    mid_pt = (length - 1) / 2
    assert 2 * mid_pt + 1 == length, 'Funny length'
    
    start = mid_pt - flank_size
    data = data[:, start: start + (flank_size * 2 + 1)]
    return data, ctl

def get_control_counts(seq_array, chosen_base, step, flank_size, sample_indices=None, circle_range=None):
    """returns the counts array for controls, more memory efficient"""
    counts = get_zero_counts((seq_array.shape[0],4), float)
    if sample_indices is None:
        sample_indices = chosen_base_indices(seq_array, chosen_base, step)
        seq_array, sample_indices = filter_seqs_by_chosen_base(seq_array, sample_indices, 1)
    
    if circle_range is None:
        circle_range = MakeCircleRange(seq_array.shape[1], flank_size)
    
    for i, v in enumerate(sample_indices):
        r = choice(v)
        indices = circle_range(r)
        counts.add_seq(seq_array[i].take(indices))
    
    return counts
