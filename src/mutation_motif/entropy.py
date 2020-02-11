"""measures of entropy for DNA sequences that are represented as numpy arrays.

Assumes bases recoded to ints in range 0 <= b < 4"""
from warnings import filterwarnings

from numpy import array, errstate, isnan, log2

from mutation_motif.util import is_valid

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-2020, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "BSD-3"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

filterwarnings("ignore", "invalid value encountered in multiply")


def get_ret(ref, ctl, pseudocount=1, check_valid=False):
    """returns relative entropy terms for axis=0

    Arguments:
        - pseudocount: amount to add to all base counts to ensure every state
          is observed
    """
    p = as_freq_matrix(ref, pseudocount=pseudocount, check_valid=check_valid)
    q = as_freq_matrix(ctl, pseudocount=pseudocount, check_valid=check_valid)
    # relative entropy terms
    with errstate(divide="ignore"):
        ret = p * log2(p / q)
    return ret


def as_freq_matrix(data, pseudocount=0, check_valid=False):
    """returns data as a frequency matrix, using pseudocounts to adjust"""
    if check_valid:
        assert is_valid(data)

    total = data.shape[0]
    base_counts = []
    for i in range(4):
        counts = (data == i).sum(axis=0) + pseudocount
        base_counts.append(counts)

    base_counts = array(base_counts, dtype=float)

    # as frequencies
    p = base_counts / (total + pseudocount * 4)
    return p


def get_entropy_terms(data, pseudocount=0, check_valid=False, freq_matrix=False):
    """returns Shannons entropy terms for axis=0"""
    if freq_matrix:
        p = data
    else:
        p = as_freq_matrix(data, pseudocount=pseudocount, check_valid=check_valid)

    with errstate(divide="ignore"):
        et = -p * log2(p)
    return et


def get_mit(data, pseudocount=0, check_valid=False, freq_matrix=False):
    """returns MI terms for axis=0"""
    entropy_terms = get_entropy_terms(
        data, pseudocount=pseudocount, check_valid=check_valid, freq_matrix=freq_matrix
    )
    mit = 0.5 - entropy_terms
    mit[isnan(mit)] = 0
    return mit
