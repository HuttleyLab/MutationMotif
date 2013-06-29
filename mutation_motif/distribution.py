"""make distribution of RE"""
from __future__ import division
import numpy
from mutation_motif import profile, entropy

def make_distribution(data, chosen_base, step, flank_size, num_reps=100):
    """returns maximum RE per position from randomised data"""
    midpoint = (data.shape[1] - 1) // 2
    assert midpoint * 2 + 1 == data.shape[1]
    randomised_data = profile.get_control(data, chosen_base, step, midpoint)
    
    start = midpoint - flank_size
    ref = randomised_data[:, start: start + (flank_size * 2 + 1)]
    
    # control generator
    make_control = profile.MakeControl(randomised_data, chosen_base, step,
                                            flank_size)
    records = []
    for i in range(num_reps):
        ctl = make_control()
        rets = entropy.get_ret(ref, ctl)
        res = rets.sum(axis=0)
        records.append(res)
    
    records = numpy.array(records, dtype=float)
    return records.max(axis=0)
