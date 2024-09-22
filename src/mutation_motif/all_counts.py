"""combines counts from each mutation direction into a single table"""

import re
from collections import Counter

import click

_directions = [
    "AtoC",
    "AtoG",
    "AtoT",
    "CtoA",
    "CtoG",
    "CtoT",
    "GtoA",
    "GtoC",
    "GtoT",
    "TtoA",
    "TtoC",
    "TtoG",
]
direction = re.compile(f'({"|".join(_directions)})')


def check_found_filenames(filenames):
    """check the number of filenames and that they include the direction"""
    found = Counter()
    for fn in filenames:
        d = direction.findall(fn)
        found.update(d)

    total = sum(found.values())
    if total != 12 or set(found) != set(_directions):
        msg = (
            f"ERROR: counts_pattern did not identify 12 files -- {filenames}\n"
            "Note that each file must contain a single direction pattern"
            ", e.g. CtoT, AtoG",
        )
        click.secho(msg, fg="red")
        exit(-1)
