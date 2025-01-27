from importlib import resources

import pytest
from cogent3 import DNA, load_aligned_seqs, make_table
from numpy import array
from numpy.testing import assert_array_equal

from mutation_motif.util import (
    array_to_str,
    get_grid_config,
    just_nucs,
    make_consistent_direction_style,
    seqs_to_array,
)


@pytest.fixture(scope="session")
def aln(DATA_DIR):
    return load_aligned_seqs(
        DATA_DIR / "just_nuc.fasta",
        array_align=True,
        moltype=DNA,
    )


def test_just_nucs(aln):
    """in the input alignment profile,
    seq0, seq2 to seq5 contain no N/-,
        should all pass the test.
    seq6 contains 2 Ns, seq7 contains 2 gaps,
    seq1 and seq8 contains both N and gap,
    seq1, seq6 to seq8 are expected to be eliminated by the code
    """
    nucs = just_nucs(aln.array_seqs)

    assert_array_equal(
        nucs,
        [
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0],
        ],
    )


@pytest.fixture(scope="session")
def d_aln(DATA_DIR):
    return load_aligned_seqs(
        DATA_DIR / "load_seqs_to_array.fasta",
        array_align=True,
        moltype=DNA,
    )


def test_seqs_to_array(d_aln):
    """in the input alignment profile,
    seq0, seq2 to seq5 contain no N/-,
        should all pass the test.
    seq6 contains 2 Ns, seq7 contains 2 gaps,
    seq1 and seq8 contains both N and gap,

    expect to convert seq0, seq2 to seq5 into numpy array
    """
    data = seqs_to_array(d_aln)
    expect = array(
        [
            (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
        ],
    )

    assert_array_equal(data, expect)


@pytest.fixture(scope="session")
def seq_array():
    return array([(2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 0)])


def test_array_to_str(seq_array):
    """convert numpy array back to DNA sequence"""
    dna_str = array_to_str(seq_array)
    expect = ["ATCAACATATAAAAAGGAAAT"]

    assert dna_str == expect


@pytest.fixture(scope="session")
def grid_cfg_path():
    return resources.files("mutation_motif") / "cfgs/grid.cfg"


def test_grid_cfg_consistency(tmp_path, grid_cfg_path):
    """fails if num rows/cols don't match paths sections"""
    cfg = grid_cfg_path.read_text()
    out = tmp_path / "grid.cfg"
    out.write_text(cfg.replace("num_cols=2", "num_cols=1"))
    with pytest.raises(ValueError):
        get_grid_config(str(out))


def test_grid_cfg(tmp_path, grid_cfg_path):
    """exercising parser"""
    resources.files("mutation_motif") / "cfgs/grid.cfg"
    cfg = grid_cfg_path.read_text()
    out = tmp_path / "grid.cfg"
    out.write_text(cfg)
    cfg = get_grid_config(str(out))


@pytest.fixture(
    scope="session",
    params=(["AtoC", "AtoG", "TtoC", "AtoT"], ["A>C", "A>G", "T>C", "A>T"]),
)
def table(request):
    data = {
        "count": [1599, 1153, 896, 711],
        "direction": request.param,
    }
    return make_table(data=data)


def test_make_consistent_style(table):
    """correctly converts X>Y to XtoY"""
    result = make_consistent_direction_style(table)
    assert result.columns["direction"].tolist() == ["AtoC", "AtoG", "TtoC", "AtoT"]


@pytest.fixture(
    scope="session",
    params=(["A:C", "A:C", "A:C", "A:C"], ["A to C", "A to C", "A to C", "A to C"]),
)
def bad_table(request):
    data = {
        "count": [1599, 1153, 896, 711],
        "direction": request.param,
    }
    return make_table(data=data)


def test_make_consistent_style_errors(bad_table):
    with pytest.raises((ValueError, AssertionError)):
        _ = make_consistent_direction_style(bad_table)
