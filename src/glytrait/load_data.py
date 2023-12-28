import csv
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import NewType

import pandas as pd

from glytrait.exception import FileTypeError, FileFormatError, NotEnoughGroupsError
from glytrait.glycan import (
    parse_structures,
    StructureDict,
    parse_compositions,
    CompositionDict,
)

AbundanceTable = NewType("AbundanceTable", pd.DataFrame)
"""Abundance table type.
The index are glycans and the columns are samples.
The abundance table could only be returned by `load_abundance_table` function.
"""

GroupSeries = NewType("GroupSeries", pd.Series)
"""Group series type.
The index are the sample names and the values are groups.
"""


def load_abundance_table(filepath: str) -> AbundanceTable:
    """Load abundance table from filepath.

    Args:
        filepath: Path to abundance table file.

    Returns:
        AbundanceTable: Abundance table.

    Raises:
        FileTypeError: If the file is not a csv file.
        FileNotFoundError: If the file does not exist.
        FileFormatError: If the file format is incorrect.
            This includes: (1) the file is not numeric; (2) the file has duplicated
            index (glycans); (3) the file has duplicated columns (samples); (4) the
            file does not have a "GlycanID" column.
    """
    _check_file_type(filepath, "csv")
    _check_exist(filepath)
    _check_columns("abundance table", filepath, ["GlycanID"])
    df = pd.read_csv(filepath, index_col="GlycanID")
    _check_numeric("abundance table", df)
    _check_duplicated_index("abundance table", df)
    _check_duplicated_columns("abundance table", df)
    return AbundanceTable(df)


def load_structures(filepath: str) -> StructureDict:
    """Load structures from filepath.

    Args:
        filepath: Path to structures file.

    Returns:
        StructureDict: A dict with glycan ids with keys and
            `Structure` objects as values.

    Raises:
        FileTypeError: If the file is not a csv file.
        FileNotFoundError: If the file does not exist.
        FileFormatError: If the file format is incorrect.
            This includes: (1) the file has duplicated glycan ids; (2) the file has
            duplicated structures; (3) the file does not have a "GlycanID" column
            or a "Structure" column; (4) the file has extra columns.
        StructureParseError: If the structure cannot be parsed.
            All glycan ids not parsed will be listed in the error message.
    """
    _check_exist(filepath)
    if Path(filepath).is_dir():
        return _load_structures_from_dir(filepath)
    else:
        return _load_structures_from_csv(filepath)


def _load_structures_from_csv(filepath: str) -> StructureDict:
    """Load structures from a csv file."""
    _check_file_type(filepath, "csv")
    _check_columns("glycan structures", filepath, ["GlycanID", "Structure"], only=True)
    ids, strings = _load_glycans_from_csv(filepath, "Structure")
    _check_duplicated_items("glycan ids", ids)
    _check_duplicated_items("glycan structures", strings)
    return parse_structures(zip(ids, strings))


def _load_structures_from_dir(dirpath: str) -> StructureDict:
    """Load structures from a directory."""
    ids, strings = _load_glycans_from_dir(dirpath, "glycoct")
    # glycan ids are not checked because they are filenames
    _check_duplicated_items("glycan structures", strings)
    return parse_structures(zip(ids, strings))


def load_compositions(filepath: str) -> CompositionDict:
    """Load compositions from filepath.

    Args:
        filepath: Path to compositions file.

    Returns:
        CompositionDict: A dict with glycan ids with keys and
            `Composition` objects as values.

    Raises:
        FileTypeError: If the file is not a csv file.
        FileNotFoundError: If the file does not exist.
        FileFormatError: If the file format is incorrect.
            This includes: (1) the file has duplicated glycan ids; (2) the file has
            duplicated compositions; (3) the file does not have a "GlycanID" column
            or a "Composition" column; (4) the file has extra columns.
        CompositionParseError: If the composition cannot be parsed.
            All glycan ids not parsed will be listed in the error message.
    """
    _check_file_type(filepath, "csv")
    _check_exist(filepath)
    _check_columns(
        "glycan compositions", filepath, ["GlycanID", "Composition"], only=True
    )
    ids, strings = _load_glycans_from_csv(filepath, "Composition")
    _check_duplicated_items("glycan ids", ids)
    _check_duplicated_items("glycan compositions", strings)
    return parse_compositions(zip(ids, strings))


def load_groups(filepath: str) -> GroupSeries:
    """Load groups from filepath.

    Args:
        filepath: Path to groups file.

    Returns:
        GroupSeries: A pandas Series with sample names as index and groups as values.

    Raises:
        FileTypeError: If the file is not a csv file.
        FileNotFoundError: If the file does not exist.
        FileFormatError: If the file format is incorrect.
            This includes: (1) the file does not have a "Group" column or a
            "Sample" column; (2) the file has extra columns; (3) the file has
            duplicated samples.
    """
    _check_file_type(filepath, "csv")
    _check_exist(filepath)
    _check_columns("groups", filepath, ["Group", "Sample"], only=True)
    s = pd.read_csv(filepath, index_col="Sample").squeeze()
    _check_duplicated_items("samples", s.index)
    _check_group_number(s)
    return GroupSeries(s)


def _check_file_type(filepath: str, file_type: str) -> None:
    if Path(filepath).suffix != f".{file_type}":
        raise FileTypeError(f"File {filepath} is not a {file_type} file.")


def _check_exist(filepath: str) -> None:
    if not Path(filepath).exists():
        raise FileNotFoundError(f"File {filepath} does not exist.")


def _check_numeric(df_name: str, df: pd.DataFrame) -> None:
    if not df.map(lambda x: isinstance(x, (int, float))).all().all():
        raise FileFormatError(f"The {df_name} should be numeric.")


def _check_duplicated_index(df_name: str, df: pd.DataFrame) -> None:
    if df.index.duplicated().any():
        raise FileFormatError(f"The {df_name} should not have duplicated index.")


def _check_duplicated_columns(df_name: str, df: pd.DataFrame) -> None:
    for col in df.columns:
        if col.endswith(".1"):
            raise FileFormatError(f"The {df_name} should not have duplicated columns.")


def _check_duplicated_items(item_name: str, to_check: Sequence[str]) -> None:
    if len(set(to_check)) != len(to_check):
        raise FileFormatError(f"The {item_name} should not be duplicated.")


def _check_columns(
    filename: str, filepath: str, columns: Iterable[str], only: bool = False
) -> None:
    def _format_columns() -> str:
        return ", ".join(f"'{col}'" for col in columns)

    with open(filepath) as f:
        header = f.readline().strip().split(",")

    if only:
        if not set(header) == set(columns):
            msg = f"The {filename} should only have columns {_format_columns()}."
            raise FileFormatError(msg)
    else:
        if not set(columns) <= set(header):
            msg = f"The {filename} should have columns {_format_columns()}."
            raise FileFormatError(msg)


def _load_glycans_from_csv(
    filepath: str, string_col: str
) -> tuple[list[str], list[str]]:
    """Helper function for `load_structures` and `load_compositions`."""
    ids: list[str] = []
    strings: list[str] = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.append(row["GlycanID"])
            strings.append(row[string_col])
    return ids, strings


def _load_glycans_from_dir(dirpath: str, suffix: str) -> tuple[list[str], list[str]]:
    """Helper function for `load_structures`."""
    ids: list[str] = []
    strings: list[str] = []
    for filepath in Path(dirpath).glob(f"*.{suffix}"):
        ids.append(filepath.stem)
        strings.append(filepath.read_text())
    return ids, strings


def _check_group_number(group_series: Iterable[str]) -> None:
    if len(set(group_series)) < 2:
        raise NotEnoughGroupsError("There should be at least two groups.")
