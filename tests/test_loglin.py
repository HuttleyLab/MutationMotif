import numpy
import pytest
from cogent3 import make_table

from mutation_motif.log_lin import position_effect, spectra_difference
from mutation_motif.util import load_table_from_delimited_file, spectra_table


@pytest.fixture(scope="session")
def single_counts():
    rows = [
        ["M", "A", 250],
        ["M", "C", 247],
        ["M", "G", 200],
        ["M", "T", 310],
        ["R", "A", 268],
        ["R", "C", 227],
        ["R", "G", 204],
        ["R", "T", 308],
    ]

    return make_table(header=("mut", "base", "count"), rows=rows)


@pytest.fixture(scope="session")
def spectra(DATA_DIR):
    path = DATA_DIR / "counts-combined.txt"
    table = load_table_from_delimited_file(path)
    return spectra_table(table, "strand")


def test_position_effect(single_counts):
    expected_deviance = 1.515818109583372  # from R
    r = position_effect(single_counts)
    numpy.testing.assert_almost_equal(r.deviance, expected_deviance, decimal=10)
    numpy.testing.assert_almost_equal(r.pvalue, 0.678624360576904)
    assert r.nfp == 3
    # ret terms computed using rpy2
    rpy = make_table(
        header=("mut", "base", "ret"),
        data=[
            ("M", "A", -0.004390161846736536),
            ("M", "C", 0.0050685522774148875),
            ("M", "G", -0.0009881162714170135),
            ("M", "T", 0.0004973269018114649),
            ("R", "A", 0.0045454767832316615),
            ("R", "C", -0.004858986131784593),
            ("R", "G", 0.0009979485970079967),
            ("R", "T", -0.0004957200241533225),
        ],
    )
    rpy = rpy.sorted(columns=["mut", "base"])
    # computed using statsmodels
    mm = make_table(data_frame=r.df)
    mm = mm.sorted(columns=["mut", "base"])
    numpy.testing.assert_almost_equal(mm.columns["ret"], rpy.columns["ret"], decimal=10)


def test_spectra_analysis(spectra):
    spectra = make_table(
        data={
            "count": [190, 190, 760, 763, 190, 190],
            "direction": ["AtoT", "AtoT", "AtoG", "AtoG", "AtoC", "AtoC"],
            "strand": ["+", "-", "+", "-", "+", "-"],
        },
    )
    expected_deviance = 0.001967210709999
    expected_pvale = 0.999016878226183
    expected_df = 2
    r = spectra_difference(spectra, "strand")
    numpy.testing.assert_almost_equal(r.deviance, expected_deviance, decimal=10)
    numpy.testing.assert_almost_equal(r.pvalue, expected_pvale, decimal=10)
    assert r.nfp == expected_df
