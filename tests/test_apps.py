import os
import shutil
from unittest import TestCase, main

from click.testing import CliRunner

from cogent3 import load_table
from mutation_motif.all_counts import main as all_count_main
from mutation_motif.aln_to_counts import main as aln_to_counts_main
from mutation_motif.draw import main as draw_main
from mutation_motif.mutation_analysis import main as mut_main
from mutation_motif.util import makedirs


class TestCounting(TestCase):
    dirname = "_delme_counts"

    def setUp(self) -> None:
        shutil.rmtree(self.dirname, ignore_errors=True)

    def test_all_counts(self):
        """exercising all_acounts"""
        runner = CliRunner()
        # should fail, as data files not in this directory
        r = runner.invoke(all_count_main, ["-cdata/*.txt", "-o%s" % self.dirname])
        self.assertNotEqual(r.exit_code, 0)
        r = runner.invoke(
            all_count_main, ["-cdata/directions/*.txt", "-o%s" % self.dirname]
        )
        # should produce directory containing two files
        dirlist = os.listdir(self.dirname)
        self.assertEqual(
            set(dirlist), set(["combined_counts.txt", "combined_counts.log"])
        )
        # check the contents of combined_counts
        counts = load_table(os.path.join(self.dirname, "combined_counts.txt"), sep="\t")
        # 4**4 nbrs x 12 mutations x 2 (M/R groups) = 6144
        counts = load_table(os.path.join(self.dirname, "combined_counts.txt"), sep="\t")
        self.assertEqual(counts.shape[0], 6144)
        shutil.rmtree(self.dirname)

    def test_aln_to_counts(self):
        """exercising aln_to_counts"""
        makedirs(self.dirname)
        runner = CliRunner()
        # should fail, as data files not in this directory
        r = runner.invoke(
            aln_to_counts_main,
            [
                "-adata/sample_AtoC.fasta",
                "-o%s" % self.dirname,
                "-f1",
                "--direction=AtoC",
                "-S111",
                "-F",
            ],
        )
        dirlist = os.listdir(self.dirname)
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(set(dirlist), set(["sample_AtoC.txt", "sample_AtoC.log"]))
        counts = load_table(os.path.join(self.dirname, "sample_AtoC.txt"), sep="\t")
        # two columns with pos, two groups giving shape=2*16
        self.assertEqual(counts.shape[0], 32)
        shutil.rmtree(self.dirname)


class TestMutationAnalysis(TestCase):
    dirname = "_delme"

    def tearDown(self) -> None:
        shutil.rmtree(self.dirname, ignore_errors=True)

    def test_nbr(self):
        """exercising, making sure output generated"""
        runner = CliRunner()
        r = runner.invoke(
            mut_main, ["nbr", "-1data/counts-CtoT.txt", "-o%s" % self.dirname]
        )
        self.assertEqual(r.exit_code, 0, r.exception)
        # expect the following file names
        fnames = [
            "1.json",
            "1.pdf",
            "2.json",
            "2.pdf",
            "3.json",
            "3.pdf",
            "4.json",
            "4.pdf",
            "summary.txt",
            "summary.pdf",
            "analysis.log",
        ]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)

    def test_nbr_ssym(self):
        """exercising, nbr strand symmetric analysis"""
        runner = CliRunner()
        r = runner.invoke(
            mut_main,
            [
                "nbr",
                "-1data/counts-CtoT-ss.txt",
                "-o%s" % self.dirname,
                "--strand_symmetry",
            ],
        )
        self.assertEqual(r.exit_code, 0)
        # expect the following file names
        fnames = [
            "1.json",
            "1.pdf",
            "2.json",
            "2.pdf",
            "3.json",
            "3.pdf",
            "4.json",
            "4.pdf",
            "summary.txt",
            "summary.pdf",
            "analysis.log",
        ]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)

    def test_spectra(self):
        """exercising spectra analysis code"""
        runner = CliRunner()
        r = runner.invoke(
            mut_main,
            [
                "spectra",
                "-1data/auto_intergen_combined_counts.txt",
                "-2data/auto_intron_combined_counts.txt",
                "-o%s" % self.dirname,
            ],
        )
        self.assertEqual(r.exit_code, 0)

        # expect the following file names
        fnames = [
            "spectra_analysis.json",
            "spectra_analysis.log",
            "spectra_summary.txt",
        ]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)

    def test_spectra_ssym(self):
        """exercising spectra analysis code with strand symmetry"""
        runner = CliRunner()
        r = runner.invoke(
            mut_main,
            [
                "spectra",
                "-1data/counts-combined.txt",
                "-o%s" % self.dirname,
                "--strand_symmetry",
            ],
        )
        self.assertEqual(r.exit_code, 0)

        # expect the following file names
        fnames = [
            "spectra_analysis.json",
            "spectra_analysis.log",
            "spectra_summary.txt",
        ]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)


class TestDrawGrid(TestCase):
    dirname = "_delme"

    def tearDown(self) -> None:
        shutil.rmtree(self.dirname, ignore_errors=True)

    def test_spectra_grid(self):
        """exercising draw spectra grid"""
        # first
        runner = CliRunner()
        r = runner.invoke(
            draw_main,
            [
                "spectra-grid",
                "--figpath=%s/spectra_grid.pdf" % self.dirname,
                "--json_path=data/spectra_analysis.json",
                "--group_label=strand",
            ],
        )

        self.assertEqual(r.exit_code, 0)
        fnames = ["spectra_grid.pdf", "spectra_grid.log"]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)

    def test_grid(self):
        """exercise drawing arbitrary grid"""
        runner = CliRunner()
        r = runner.invoke(
            draw_main,
            [
                "grid",
                "--figpath=%s/grid.pdf" % self.dirname,
                "--fig_config=data/arbitrary_grid.cfg",
            ],
        )

        self.assertEqual(r.exit_code, 0)
        fnames = ["grid.pdf", "grid.log"]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)


if __name__ == "__main__":
    main()
