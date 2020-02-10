import numpy
import pandas as pd
from rpy2.robjects import Formula
from rpy2.robjects import r as R
from rpy2.robjects.vectors import DataFrame, FactorVector, IntVector, StrVector

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-2020, Gavin Huttley, Yicheng Zhu"
__credits__ = ["Gavin Huttley", "Yicheng Zhu"]
__license__ = "BSD-3"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


def as_dataframe(table):
    """returns a DataFrame instance. Requires counts to be
    [[col1, col2, col3, ..]]"""
    data = dict(list(zip(table.header, list(zip(*table.tolist())))))
    for column in data:
        if type(data[column][0]) in (str, str):
            klass = StrVector
        else:
            klass = IntVector

        data[column] = klass(data[column])

    return DataFrame(data)


def get_rpy2_factorvector_index_map(fac):
    """returns dict mapping factor vector indices to factors"""
    mapping = dict([(i + 1, f) for i, f in enumerate(fac.levels)])
    return mapping


def convert_rdf_to_pandasdf(r_df):
    """converts an rpy2 dataframe to a pandas dataframe"""
    converted = {}
    for col_name, col_vector in list(r_df.items()):
        if type(col_vector) == FactorVector:
            index_factor_map = get_rpy2_factorvector_index_map(col_vector)
            col_vector = [index_factor_map[val] for val in col_vector]
        converted[col_name] = list(col_vector)
    return pd.DataFrame(converted)


def DevianceToRelativeEntropy(N):
    """converts deviance to Relative Entropy"""
    denom = 2 * N

    def call(val):
        return val / denom

    return call


def CalcRet(dev_to_re, epsilon=1e-9):
    """factory function for computing residual relative entropy terms.

    dev_to_re is a function for converting a deviance to relative entropy"""

    def call(obs, exp):
        result = []
        for i in range(len(obs)):
            o, e = obs[i], exp[i]
            e = e or epsilon  # to avoide zero division
            o = o or epsilon  # avoid zero counts
            ret = dev_to_re(2 * o * numpy.log(o / e))
            result.append(ret)
        return result

    return call


def position_effect(counts_table, group_label=None, test=False):
    """returns total relative entropy, degrees of freedom and stats

    fit's a log-lin model that excludes only the full interaction term

    Arguments:
        - group_label: name of column containing group data
    """
    num_pos = sum(1 for c in counts_table.header if c.startswith("base"))
    assert 1 <= num_pos <= 4, "Can only handle 4 positions"

    if num_pos == 1:
        columns = ["mut", "base", "count"]
    else:
        columns = ["mut"] + ["base%d" % (i + 1) for i in range(num_pos)] + ["count"]

    # handle groups
    if group_label and group_label in counts_table.header:
        columns.insert(0, group_label)

    factors = columns[:-1]
    formula = " - ".join([" * ".join(factors), " : ".join(factors)])
    formula = "count ~ %s" % formula
    null = Formula(formula)
    if test:
        print(formula)

    counts_table = counts_table.get_columns(columns)
    d = as_dataframe(counts_table)

    f = R.glm(null, data=d, family="poisson")
    f_attr = dict(list(f.items()))
    dev = f_attr["deviance"][0]
    df = f_attr["df.residual"][0]

    collated = convert_rdf_to_pandasdf(f_attr["data"])
    collated["fitted"] = list(f_attr["fitted.values"])
    dev_to_re = DevianceToRelativeEntropy(collated["count"].sum())
    calc_ret = CalcRet(dev_to_re)
    total_re = dev_to_re(dev)

    collated["ret"] = calc_ret(collated["count"], collated["fitted"])
    collated = collated.reindex(columns + ["fitted", "ret"], axis=1)
    collated = collated.sort_values(by=columns[:-1])
    return total_re, dev, df, collated, formula


def spectra_difference(counts_table, group_label, test=False):
    """group_label is the column name for category"""
    # we compare direction between group
    columns = ["count", "direction", group_label]
    assert set(columns) <= set(counts_table.header)
    formula = "count ~ direction + %s" % group_label
    null = Formula(formula)
    if test:
        print(formula)

    counts_table = counts_table.get_columns(columns)
    d = as_dataframe(counts_table)
    f = R.glm(null, data=d, family="poisson")
    f_attr = dict(list(f.items()))
    dev = f_attr["deviance"][0]
    df = f_attr["df.residual"][0]

    collated = convert_rdf_to_pandasdf(f_attr["data"])
    collated["fitted"] = list(f_attr["fitted.values"])
    dev_to_re = DevianceToRelativeEntropy(collated["count"].sum())
    calc_ret = CalcRet(dev_to_re)
    total_re = dev_to_re(dev)

    collated["ret"] = calc_ret(collated["count"], collated["fitted"])
    collated = collated.reindex(columns + ["fitted", "ret"], axis=1)
    collated = collated.sort_values(by=columns[:-1])
    return total_re, dev, df, collated, formula
