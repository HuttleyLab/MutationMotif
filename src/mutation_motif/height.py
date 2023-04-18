from warnings import filterwarnings

from numpy import errstate, fabs, isnan


filterwarnings("ignore", "invalid value encountered.*")


def get_mi_char_heights(freq_matrix, mi):
    """returns char height for each position in aln."""
    result = freq_matrix * mi
    result[isnan(result)] = 0
    return result


def get_re_char_heights(rets, re_positionwise=None):
    """returns character proportions for each position, can be negative

    Arguments:
        - rets: position wise relative entropy terms
        - re_positionwise: total relative entropy for each position. If None,
          computed as the sum of the position wise rets. NOTE: for using ret's
          from a log-lin model, provide these values
    """
    pwise_intervals = fabs(rets).sum(axis=0)  # span of column terms
    with errstate(divide="ignore"):
        normalised = rets / pwise_intervals

    if re_positionwise is None:
        re_positionwise = rets.sum(axis=0)

    result = normalised * re_positionwise
    result[isnan(result)] = 0
    return result
