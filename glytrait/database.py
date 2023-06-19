from __future__ import annotations

from importlib.resources import as_file, files
from typing import Literal, Iterable

from glytrait.glycan import NGlycan
from glytrait.io import read_structure

built_in_db = {
    "serum": files("glytrait.resources").joinpath("serum_structures.csv"),
    "IgG": files("glytrait.resources").joinpath("IgG_structures.csv"),
}


def load_default(
    db: Literal["serum", "IgG"], compositions: Iterable[str]
) -> list[NGlycan]:
    """Return a list of default structures.

    Args:
        db: The database to use. Either "serum" or "IgG".
        compositions: The glycan compositions to search for.

    Returns:
        A list of default structures, in the order of the compositions.
    """
    db_file = built_in_db[db]
    with as_file(db_file) as db_filename:
        return read_structure(str(db_filename), compositions)
