import os
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner
from cogent3 import load_table

from mutation_motif.cli import main as mut_main

RUNNER = CliRunner()


def test_all_counts_fail(tmp_dir, DATA_DIR):
    # should fail, as data files not in this directory
    data_glob = DATA_DIR / "*.txt"
    r = RUNNER.invoke(
        mut_main,
        ["prep-spectra", f"-c{data_glob}", f"-o{tmp_dir!s}"],
        catch_exceptions=False,
    )
    assert r.exit_code != 0


def test_all_counts(tmp_dir, DATA_DIR):
    path = DATA_DIR / "directions/*.txt"
    r = RUNNER.invoke(
        mut_main,
        ["prep-spectra", f"-c{path}", f"-o{tmp_dir!s}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0
    # should produce directory containing two files
    dirlist = os.listdir(tmp_dir)
    assert set(dirlist) == {"combined_counts.txt", "combined_counts.log"}

    # check the contents of combined_counts
    counts = load_table(tmp_dir / "combined_counts.txt", sep="\t")
    # 4**4 nbrs x 12 mutations x 2 (M/R groups) = 6144
    assert counts.shape[0] == 6144


def test_all_counts_1(tmp_dir, DATA_DIR):
    """exercising all_counts with strand symmetric"""
    path = DATA_DIR / "directions/*.txt"
    r = RUNNER.invoke(
        mut_main,
        ["prep-spectra", f"-c{path}", f"-o{tmp_dir}", "-s"],
    )

    # should produce directory containing two files
    assert r.exit_code == 0, r.output
    dirlist = {p.name for p in tmp_dir.glob("*")}
    assert dirlist == {"combined_counts.txt", "combined_counts.log"}

    counts = load_table(tmp_dir / "combined_counts.txt", sep="\t")
    assert "strand" in counts.header


def test_all_counts_splitdir(tmp_dir, DATA_DIR):
    path = DATA_DIR / "directions/*.txt"

    splitdir = tmp_dir / "splitdir"
    r = RUNNER.invoke(
        mut_main,
        [
            "prep-spectra",
            f"-c{path}",
            f"-o{tmp_dir}",
            "-s",
            f"--split_dir={splitdir}",
        ],
    )
    # should produce directory containing 6 txt files
    assert r.exit_code == 0, r.output
    dirlist = {p.name for p in splitdir.glob("*")}
    assert len(dirlist) == 6
    expected_size = 4**4 * 2 * 2
    tables = [load_table(p, sep="\t") for p in splitdir.glob("*")]
    assert all(t.shape[0] == expected_size for t in tables)
    assert all("strand" in t.header for t in tables)


def test_aln_to_counts(tmp_path, DATA_DIR):
    """exercising aln_to_counts"""
    path = DATA_DIR / "sample_AtoC.fasta"

    # should fail, as data files not in this directory
    r = RUNNER.invoke(
        mut_main,
        [
            "prep-nbr",
            f"-a{path}",
            f"-o{tmp_path}",
            "-f1",
            "--direction=AtoC",
            "-S111",
            "-F",
        ],
    )
    assert r.exit_code == 0, r.output
    dirlist = list(tmp_path.glob("*"))

    assert {p.name for p in dirlist} == {"sample_AtoC.txt", "sample_AtoC.log"}
    counts = load_table(tmp_path / "sample_AtoC.txt", sep="\t")
    # two columns with pos, two groups giving shape=2*16
    assert counts.shape[0] == 32


def test_nbr(tmp_path, DATA_DIR):
    """exercising, making sure output generated"""
    path = DATA_DIR / "counts-CtoT.txt"
    r = RUNNER.invoke(
        mut_main,
        ["ll-nbr", f"-1{path}", f"-o{tmp_path}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
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
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_nbr_ssym(tmp_path, DATA_DIR):
    """exercising, nbr strand symmetric analysis"""
    path = DATA_DIR / "counts-CtoT-ss.txt"
    r = RUNNER.invoke(
        mut_main,
        [
            "ll-nbr",
            f"-1{path}",
            f"-o{tmp_path}",
            "--strand_symmetry",
        ],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
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
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_spectra(tmp_path, DATA_DIR):
    """exercising spectra analysis code"""
    path1 = DATA_DIR / "auto_intergen_combined_counts.txt"
    path2 = DATA_DIR / "auto_intron_combined_counts.txt"
    r = RUNNER.invoke(
        mut_main,
        [
            "ll-spectra",
            f"-1{path1}",
            f"-2{path2}",
            f"-o{tmp_path}",
        ],
        catch_exceptions=False,
    )

    assert r.exit_code == 0, r.output

    # expect the following file names
    fnames = [
        "spectra_analysis.json",
        "spectra_analysis.log",
        "spectra_summary.txt",
    ]
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_spectra_ssym(tmp_path, DATA_DIR):
    """exercising spectra analysis code with strand symmetry"""
    path = DATA_DIR / "counts-combined.txt"
    r = RUNNER.invoke(
        mut_main,
        [
            "ll-spectra",
            f"-1{path}",
            f"-o{tmp_path}",
            "--strand_symmetry",
        ],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output

    # expect the following file names
    fnames = [
        "spectra_analysis.json",
        "spectra_analysis.log",
        "spectra_summary.txt",
    ]
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_spectra_grid(tmp_path, DATA_DIR):
    """exercising draw spectra grid"""
    path = DATA_DIR / "spectra_analysis.json"
    # first
    r = RUNNER.invoke(
        mut_main,
        [
            "draw-spectra-grid",
            f"--figpath={tmp_path}/spectra_grid.pdf",
            f"--json_path={path}",
            "--group_label=strand",
        ],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    fnames = ["spectra_grid.pdf", "spectra_grid.log"]
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_grid(tmp_path, DATA_DIR):
    """exercise drawing arbitrary grid"""
    path = DATA_DIR / "arbitrary_grid.cfg"
    r = RUNNER.invoke(
        mut_main,
        [
            "draw-grid",
            f"--figpath={tmp_path}/grid.pdf",
            f"--fig_config={path}",
        ],
        catch_exceptions=False,
    )

    assert r.exit_code == 0, r.output
    fnames = ["grid.pdf", "grid.log"]
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_nbr_app(DATA_DIR, tmp_path):
    """cl produces plots for 1-way up to 4-way plus summary"""

    data_path = tmp_path / "CtoT"
    shutil.copytree(DATA_DIR / "CtoT", data_path)
    cfg_path = Path(tmp_path) / "nbr_paths.cfg"
    cfg_path.write_text((DATA_DIR / "nbr_paths.cfg").read_text())
    r = RUNNER.invoke(
        mut_main,
        ["draw-nbr", f"-p{cfg_path}"],
    )
    assert r.exit_code == 0, r.output
    fnames = [f"{n}.pdf" for n in ("one", "two", "three", "four", "summary")]
    assert all(f.stat().st_size > 0 for f in (tmp_path / fn for fn in fnames))


def test_nbr_matrix_app(DATA_DIR, tmp_path):
    """cl produces matrix of 1-way plots"""

    data_path = tmp_path / "directions"
    shutil.copytree(DATA_DIR / "directions", data_path)
    cfg_path = Path(tmp_path) / "nbr_matrix_paths.cfg"
    shutil.copy(DATA_DIR / "nbr_matrix_paths.cfg", cfg_path)
    figpath = Path(tmp_path) / "nbr_matrix.pdf"
    r = RUNNER.invoke(
        mut_main,
        ["draw-nbr-matrix", f"--paths_cfg={cfg_path}", f"--figpath={figpath}"],
    )
    assert r.exit_code == 0, r.output
    assert figpath.exists()
    assert figpath.stat().st_size > 0


@pytest.mark.parametrize("use_freq", [False, True])
def test_mi_app(DATA_DIR, tmp_path, use_freq):
    """cl produces 1-way plot using MI"""
    data_path = DATA_DIR / "directions" / "CtoT.json"
    figpath = tmp_path / "mi.pdf"
    args = ["draw-mi", f"--json_path={data_path}", f"--figpath={figpath}"]
    args += ["--use_freq"] if use_freq else []
    r = RUNNER.invoke(
        mut_main,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert figpath.exists()
    assert figpath.stat().st_size > 0


def test_export_cfg_app(tmp_path):
    """exports sample cfg files"""

    path = tmp_path / "cfgs"
    r = RUNNER.invoke(
        mut_main,
        ["draw-export-cfg", str(path)],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    file_stats = [
        f.stat().st_size > 0 for f in (tmp_path / fn for fn in path.glob("*.cfg"))
    ]
    assert all(file_stats)
    assert file_stats
