import dataclasses
import typing

import numpy
import pandas
import statsmodels.api as sm
import statsmodels.formula.api as smf
from cogent3.util.table import Table
from scipy.stats import chi2

_poisson = sm.families.Poisson()
_glm = smf.glm


def table_to_pandas_with_factors(
    table: Table,
    factor_columns: list[str],
) -> pandas.DataFrame:
    """Returns a pandas DataFrame with factor columns converted to pandas category type"""
    df = table.to_pandas()
    for column in factor_columns:
        df[column].astype("category")
    return df


class DevianceToRelativeEntropy:
    """converts a deviance to relative entropy"""

    def __init__(self, N: int) -> None:
        self.denom = 2 * N

    def __call__(self, val: int | float) -> float:
        return val / self.denom


class CalcRet:
    """computes residual relative entropy terms."""

    def __init__(
        self,
        dev_to_re: typing.Callable[[float], float],
        epsilon: float = 1e-9,
    ):
        self.dev_to_re = dev_to_re
        self.epsilon = epsilon

    def __call__(
        self,
        obs: list[float],
        exp: list[float],
    ) -> list[float]:
        result = []
        for i in range(len(obs)):
            o, e = obs[i], exp[i]
            e = e or self.epsilon  # to avoide zero division
            o = o or self.epsilon  # avoid zero counts
            ret = self.dev_to_re(2 * o * numpy.log(o / e))
            result.append(ret)
        return result


@dataclasses.dataclass
class loglin_result:
    relative_entropy: float
    deviance: float
    nfp: int
    formula: str
    df: pandas.DataFrame
    pvalue: float


def stats_from_loglin(
    counts_table: Table,
    factors: list[str],
    formula: str,
) -> loglin_result:
    df = table_to_pandas_with_factors(counts_table, factors)
    model = _glm(formula=formula, data=df, family=_poisson).fit()
    dev = model.deviance
    nfp = model.df_resid
    dev_to_re = DevianceToRelativeEntropy(counts_table.sum_columns("count"))
    calc_ret = CalcRet(dev_to_re)
    total_re = dev_to_re(dev)
    df["fitted"] = model.predict()
    df["ret"] = calc_ret(df["count"], df["fitted"])
    return loglin_result(
        relative_entropy=float(total_re),
        deviance=float(dev),
        nfp=int(nfp),
        formula=formula,
        df=df,
        pvalue=chi2.sf(dev, nfp),
    )


def position_effect(
    counts_table: Table,
    group_label: str | None = None,
    test: bool = False,
) -> loglin_result:
    """fit's a log-linear model for a single mutation direction

    Parameters
    ----------
    counts_table
        table of counts for mutated bases and their controls
    group_label
        will group data by this column and test for a difference between the
        groups
    test
        verbose output

    Notes
    -----
    The method is assuming the counts are for a single mutation direction.
    The log-linear model excludes only the full interaction term.
    """
    num_pos = sum(bool(c.startswith("base")) for c in counts_table.header)
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
    formula = f"count ~ {formula}"
    if test:
        print(formula)

    counts_table = counts_table.get_columns(columns)
    return stats_from_loglin(counts_table, factors, formula)


def spectra_difference(
    counts_table: Table,
    group_label: str,
    test: bool = False,
) -> loglin_result:
    """fits a log-linear model for equivalence of spectra between groups

    Parameters
    ----------
    counts_table
        table of counts for mutation outcomes from one starting base, split by
        group where group can be, for example, "strand".
    group_label
        the column containing the group factor labels
    test
        verbose output

    Notes
    -----
    For a A mmutation spectra is the distribution of mutations across the four
    The log-linear model excludes only the full interaction term.
    """
    # we compare direction between group
    columns = ["count", "direction", group_label]
    assert set(columns) <= set(counts_table.header)
    formula = f"count ~ direction + {group_label}"
    if test:
        print(formula)

    counts_table = counts_table.get_columns(columns)
    return stats_from_loglin(counts_table, columns[1:], formula)
