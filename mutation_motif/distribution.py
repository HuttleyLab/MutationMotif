"""make distribution of RE"""

from __future__ import division
import numpy
from cogent.maths.stats.util import UnsafeNumbers
from mutation_motif import profile, entropy

def get_re_distribution(observed, chosen_base, step, flank_size, make_control=None, num_reps=100):
    """returns positionwise RE per replicate with RE generation done using random
    controls on each iteration"""
    # TODO move the chosen_base flank_size midpoint logic out of function.
    # Restricts utility to centred mutation motifs only
    midpoint = (observed.shape[1] - 1) // 2
    assert midpoint * 2 + 1 == observed.shape[1]
    start = midpoint - flank_size
    ref = observed[:, start: start + (flank_size * 2 + 1)]
    
    # control generator
    if make_control is None:
        make_control = profile.MakeControl(observed, chosen_base, step,
                                            flank_size)
    records = []
    for i in range(num_reps):
        ctl = make_control()
        rets = entropy.get_ret(ref, ctl)
        res = rets.sum(axis=0)
        records.append(res)
    
    records = numpy.array(records, dtype=float)
    return records

def get_max_re_distribution(data, chosen_base, step, flank_size, make_control=None, num_reps=100):
    """returns maximum RE per position from randomised data"""
    midpoint = (data.shape[1] - 1) // 2
    assert midpoint * 2 + 1 == data.shape[1]
    randomised_data = profile.get_control(data, chosen_base, step, midpoint)
    
    records = get_re_distribution(randomised_data, chosen_base, step,
                          flank_size, make_control=make_control, num_reps=100)
    return records.max(axis=0)

def get_resampled_re(obs, ctl, get_et=entropy.get_ret, num_reps=100):
    """returns per sample positionwise RE from resampling (with replacement) of obs/ctl
    
    random indices are drawn, with replacement, and used to samples from both
    the provided obs and ctl arrays"""
    def calc_stat(obs, ctl):
        def call(indices):
            work_obs = obs.take(indices, axis=0)
            work_ctl = ctl.take(indices, axis=0)
            return entropy.get_ret(work_obs, work_ctl).sum(axis=0)
        
        return call
    
    num_seqs = obs.shape[0]
    records = []
    for i in range(num_reps):
        indices = numpy.random.random_integers(0, high=num_seqs-1, size=num_seqs)
        sampled_obs = obs.take(indices, axis=0)
        sampled_ctl = ctl.take(indices, axis=0)
        rets = entropy.get_ret(sampled_obs, sampled_ctl)
        res = rets.sum(axis=0)
        records.append(res)
    
    records = numpy.array(records, dtype=float)
    return records

def estimate_re_ci(data, chosen_base, step, flank_size, make_control=None, num_reps=100, lo_quantile=0.05, hi_quantile=0.95):
    """returns (lower, upper) quantiles for RE from using get_re_distribution"""
    records = get_re_distribution(data, chosen_base, step,
                          flank_size, make_control=make_control, num_reps=100)
    uppers = []
    lowers = []
    for i in range(records.shape[1]):
        position = UnsafeNumbers(records[:,i])
        upper = position.quantile(hi_quantile)
        lower = position.quantile(lo_quantile)
        uppers.append(upper)
        lowers.append(lower)
    
    return numpy.array(lowers), numpy.array(uppers)

def bootstrapped_quantiles(obs, ctl, get_et=entropy.get_ret, num_reps=100, lo_quantile=0.05, hi_quantile=0.95):
    '''returns lower, upper confidence
    
    Arguments:
        - get_et: function that returns the entropy terms. Defaults to re terms
        - num_reps: number of bootstrap samples (with replacement)
        - lo_quantile: lower bound for confindence interval
        - hi_quantile: higher bound for confindence interval
        '''
    records = get_resampled_re(obs, ctl, get_et=get_et, num_reps=num_reps)
    uppers = []
    lowers = []
    for i in range(records.shape[1]):
        position = UnsafeNumbers(records[:,i])
        uppers.append(position.quantile(hi_quantile))
        lowers.append(position.quantile(lo_quantile))
    
    return numpy.array(lowers), numpy.array(uppers)
