from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3 import DNA, load_aligned_seqs, make_table
from numpy import array
from numpy.testing import assert_array_equal
from pkg_resources import resource_filename

from mutation_motif.util import (
    array_to_str,
    get_grid_config,
    just_nucs,
    make_consistent_direction_style,
    seqs_to_array,
)


class TestJustNucs(TestCase):
    aln = load_aligned_seqs("data/just_nuc.fasta", array_align=True, moltype=DNA)

    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,
        seq1, seq6 to seq8 are expected to be eliminated by the code
        """
        nucs = just_nucs(self.aln.array_seqs)

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


class TestAlignSnpAnnotation(TestCase):
    d_aln = load_aligned_seqs(
        "data/load_seqs_to_array.fasta", array_align=True, moltype=DNA
    )

    def test_seqs_to_array(self):
        """in the input alignment profile,
        seq0, seq2 to seq5 contain no N/-,
            should all pass the test.
        seq6 contains 2 Ns, seq7 contains 2 gaps,
        seq1 and seq8 contains both N and gap,

        expect to convert seq0, seq2 to seq5 into numpy array
        """
        data = seqs_to_array(self.d_aln)
        expect = array(
            [
                (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                (2, 0, 1, 2, 3, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 0, 2, 2, 3, 0),
            ]
        )

        assert_array_equal(data, expect)


class TestAlignSnpAnnotation2(TestCase):
    seq_array = array([(2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 0)])

    def test_array_to_str(self):
        """convert numpy array back to DNA sequence"""
        dna_str = array_to_str(self.seq_array)
        expect = ["ATCAACATATAAAAAGGAAAT"]

        self.assertEqual(dna_str, expect)


class TestCfgParsing(TestCase):
    def test_grid_cfg_consistency(self):
        """fails if num rows/cols don't match paths sections"""
        path = resource_filename("mutation_motif", f"cfgs/grid.cfg")
        cfg = Path(path).read_text()
        with TemporaryDirectory(dir=".") as dirname:
            out = Path(dirname) / "grid.cfg"
            out.write_text(cfg.replace("num_cols=2", "num_cols=1"))
            with self.assertRaises(ValueError):
                get_grid_config(str(out))

    def test_grid_cfg(self):
        """exercising parser"""
        path = resource_filename("mutation_motif", f"cfgs/grid.cfg")
        cfg = Path(path).read_text()
        with TemporaryDirectory(dir=".") as dirname:
            out = Path(dirname) / "grid.cfg"
            out.write_text(cfg)
            cfg = get_grid_config(str(out))


class TestDirectionStyle(TestCase):
    """test conversion of X>Y to AtoY"""

    def test_make_consistent_style(self):
        """correctly converts X>Y to XtoY"""
        data = {
            "count": [1599, 1153, 896, 711],
            "direction": ["AtoC", "AtoG", "TtoC", "AtoT"],
        }
        table = make_table(data=data)
        result = make_consistent_direction_style(table)
        self.assertEqual(
            result.columns["direction"].tolist(), ["AtoC", "AtoG", "TtoC", "AtoT"]
        )

        data = {
            "count": [1599, 1153, 896, 711],
            "direction": ["A>C", "A>G", "T>C", "A>T"],
        }
        table = make_table(data=data)
        result = make_consistent_direction_style(table)
        self.assertEqual(
            result.columns["direction"].tolist(), ["AtoC", "AtoG", "TtoC", "AtoT"]
        )

        # wrong separator raises
        with self.assertRaises(ValueError):
            data = {
                "count": [1599, 1153, 896, 711],
                "direction": ["A:C", "A:C", "A:C", "A:C"],
            }
            table = make_table(data=data)
            _ = make_consistent_direction_style(table)

        # wrong length for separator raises
        with self.assertRaises(AssertionError):
            data = {
                "count": [1599, 1153, 896, 711],
                "direction": ["A to C", "A to C", "A to C", "A to C"],
            }
            table = make_table(data=data)
            _ = make_consistent_direction_style(table)


if __name__ == "__main__":
    main()
