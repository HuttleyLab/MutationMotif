from numpy import array, isnan, fabs
from cogent.align.weights.util import AlnToProfile, DNA, DNA_ORDER

# this for MI
def get_mi_char_heights(freq_matrix, mi, zero_middle=True):
    """returns char height for each position in aln.
    """
    result = freq_matrix * mi
    if zero_middle:
        # zero the middle position
        mid = (result.shape[1] -1) / 2
        result[:, mid] = 0
    return result

def get_re_char_heights(rets, re_positionwise=None, zero_middle=True):
    """returns character proportions for each position, can be negative
    
    Arguments:
        - rets: position wise relative entropy terms
        - re_positionwise: total relative entropy for each position. If None, computed
          as the sum of the position wise rets. NOTE: for using ret's from
          a log-lin model, provide these values
        - zero_middle: set the middle position to 0
    """
    pwise_intervals = fabs(rets).sum(axis=0) # span of column terms
    normalised = rets / pwise_intervals
    
    if re_positionwise is None:
        re_positionwise = rets.sum(axis=0)
    
    result = normalised * re_positionwise
    result[isnan(result)] = 0
    if zero_middle:
        # zero the middle position
        mid = (rets.shape[1] - 1) / 2
        result[:,mid] = 0
    
    return result

