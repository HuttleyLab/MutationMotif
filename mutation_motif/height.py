from numpy import array, isnan, fabs
from cogent.align.weights.util import AlnToProfile, DNA, DNA_ORDER

# this for MI
def get_mi_char_heights(probs, heights, zero_middle=True):
    """returns char height for each position in aln.
    """
    result = probs * heights
    if zero_middle:
        # zero the middle position
        mid = (result.shape[1] -1 ) / 2
        result[:, mid] = 0
    return result

def get_re_char_heights(rets, zero_middle=True):
    """returns character proportions for each position, can be negative"""
    pwise_intervals = fabs(rets).sum(axis=0) # span of column terms
    normalised = rets / pwise_intervals
    re_pwise = rets.sum(axis=0)
    result = normalised * re_pwise
    result[isnan(result)] = 0
    if zero_middle:
        # zero the middle position
        mid = (rets.shape[1] -1 ) / 2
        result[:,mid] = 0
    
    return result

