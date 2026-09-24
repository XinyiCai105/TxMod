"""Locations of the input tables and of the rendered figures.

Input tables are the processed intermediate tables deposited at Zenodo
(see paper/README.md). Point TXMOD_DATA at the unpacked deposit, or place the
files in paper/data/. Files may be kept gzip-compressed; ``name.tsv`` is
resolved to ``name.tsv.gz`` when only the compressed copy is present.
Figures are written to TXMOD_OUT, default paper/output/.
"""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("TXMOD_DATA", HERE / "data"))
OUT = Path(os.environ.get("TXMOD_OUT", HERE / "output"))


def data(name, must_exist=True):
    """Path of an input table, accepting a gzip-compressed copy."""
    for p in (DATA / name, DATA / (name + ".gz")):
        if p.exists():
            return str(p)
    if must_exist:
        raise FileNotFoundError(
            f"{name} not found in {DATA}. Download the Zenodo deposit and set "
            f"TXMOD_DATA to the directory that holds its tables.")
    return str(DATA / name)


def out(name):
    """Path for a rendered figure, creating the output directory if needed."""
    OUT.mkdir(parents=True, exist_ok=True)
    return str(OUT / name)
