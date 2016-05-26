import os, shutil

import click
from click.testing import CliRunner

from cogent.util.unit_test import TestCase, main

from mutation_motif.mutation_analysis import main as mut_main


class TestMutationAnalysis(TestCase):
    dirname = "_delme"
    def test_nbr(self):
        '''exercising, making sure output generated'''
        runner = CliRunner()
        r = runner.invoke(mut_main, ["-1data/counts-CtoT.txt", "-o%s" % self.dirname, "nbr"])
        self.assertEqual(r.exit_code, 0)
        # expect the following file names
        fnames = ["1.json", "1.pdf", "2.json", "2.pdf", "3.json", "3.pdf", "4.json", "4.pdf",
                 "summary.txt", "summary.pdf", "analysis.log"]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)
        shutil.rmtree(self.dirname)

    def test_nbr_ssym(self):
        '''exercising, nbr strand symmetric analysis'''
        runner = CliRunner()
        r = runner.invoke(mut_main, ["-1data/counts-CtoT-ss.txt", "-o%s" % self.dirname,
                        "--strand_symmetry", "nbr"])
        self.assertEqual(r.exit_code, 0)
        # expect the following file names
        fnames = ["1.json", "1.pdf", "2.json", "2.pdf", "3.json", "3.pdf", "4.json", "4.pdf",
                 "summary.txt", "summary.pdf", "analysis.log"]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)
        shutil.rmtree(self.dirname)
    
    def test_spectra(self):
        """exercising spectra analysis code"""
        runner = CliRunner()
        r = runner.invoke(mut_main, ["-1data/counts-combined.txt", "-o%s" % self.dirname,
                        "--strand_symmetry", "spectra"])
        self.assertEqual(r.exit_code, 0)
        
        # expect the following file names
        fnames = ["spectra_analysis.json", "spectra_analysis.log", "spectra_summary.txt"]
        for fn in fnames:
            path = os.path.join(self.dirname, fn)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.getsize(path) > 0)
        shutil.rmtree(self.dirname)

if __name__ == '__main__':
    main()
