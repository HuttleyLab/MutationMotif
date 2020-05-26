#!/usr/bin/env python
import pathlib
import sys

from setuptools import find_packages, setup


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "0.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

# Check Python version, no point installing if unsupported version inplace
if sys.version_info < (3, 6):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError(
        "Python-3.6 or greater is required, Python-%s used." % py_version
    )

short_description = "MutationMotif"

# This ends up displayed by the installer
long_description = (
    "MutationMotif implements statistical methods for the "
    "analysis of point mutation spectra and provides sequence "
    "logo style visualisation of the fitted models. "
    "The tests evaluate the association of flanking DNA sequence "
    "with mutation direction. Statistical tests for contrasting "
    "samples are also included. The tests control for "
    "heterogeneity in sequence composition across the genome. "
    "The methods are described in "
    "`Zhu et al (2017) <https://www.ncbi.nlm.nih.gov/pubmed/27974498>`_."
)

PACKAGE_DIR = "src"

setup(
    name="mutation_motif",
    version=__version__,
    author="Gavin Huttley",
    author_email="gavin.huttley@anu.edu.au",
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/HuttleyLab/mutationmotif",
    platforms=["any"],
    license=["BSD"],
    keywords=[
        "biology",
        "genomics",
        "genetics",
        "statistics",
        "evolution",
        "bioinformatics",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "cogent3",
        "click",
        "pandas",
        "plotly",
        "psutil",
        "requests",
        "rpy2",
        "scitrack",
    ],
    # note: http://stackoverflow.com/questions/3472430/how-can-i-make-setuptools-install-a-package-thats-not-on-pypi
    # and http://stackoverflow.com/questions/17366784/setuptools-unable-to-use-link-from-dependency-links/17442663#17442663
    # changing it to http://github.com/mtai/python-gearman/tarball/master#egg=gearman-2.0.0beta instead
    entry_points={
        "console_scripts": [
            "mutation_analysis=mutation_motif.mutation_analysis:main",
            "aln_to_counts=mutation_motif.aln_to_counts:main",
            "all_counts=mutation_motif.all_counts:main",
            "mutation_draw=mutation_motif.draw:main",
        ],
    },
    packages=find_packages(where="src"),
    package_dir={"": PACKAGE_DIR},
    package_data={
        "mutation_motif": [
            "cfgs/grid.cfg",
            "cfgs/nbr.cfg",
            "cfgs/nbr_matrix.cfg",
            "cfgs/nbr_paths.cfg",
            "cfgs/spectra.cfg",
        ]
    },
)
