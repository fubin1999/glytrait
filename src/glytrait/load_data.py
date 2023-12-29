from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Optional, Literal

import pandas as pd
from attrs import define

from glytrait.data_type import (
    AbundanceTable,
    GroupSeries,
)
from glytrait.exception import FileTypeError, FileFormatError, NotEnoughGroupsError
from glytrait.glycan import (
    parse_structures,
    parse_compositions,
    StructureDict,
    CompositionDict,
    GlycanDict,
)

__all__ = [
    "GlyTraitInputData",
    "load_input_data",
    "load_abundance_table",
    "load_structures",
    "load_compositions",
    "load_groups",
]


# ===== The highest-level API =====
def load_input_data(
    *,
    abundance_file: str,
    glycan_file: str,
    mode: Literal["structure", "composition"] = "structure",
    group_file: Optional[str] = None,
) -> GlyTraitInputData:
    """Load all the input data for GlyTrait.

    Notes:
        If `group_file` is not provided, the `groups` attribute of the returned
        `GlyTraitInputData` will be `None`.

    Args:
        abundance_file: Path to the abundance table file.
        glycan_file: Path to the glycan file.
        mode: Either "structure" or "composition".
        group_file: Path to the group file. Optional.

    Returns:
        GlyTraitInputData: Input data for GlyTrait.

    Raises:
        FileTypeError: If any of the files is not a csv file.
        FileNotFoundError: If any of the files does not exist.
        FileFormatError: If any of the files has incorrect format.
            See docstrings of `load_abundance_table`, `load_structures`,
            `load_compositions`, and `load_groups` for details.
    """
    abundance_table = load_abundance_table(abundance_file)
    glycans: GlycanDict
    if mode == "structure":
        glycans = load_structures(glycan_file)
    elif mode == "composition":
        glycans = load_compositions(glycan_file)
    else:
        raise ValueError(f"Invalid mode {mode}.")
    if group_file is not None:
        groups = load_groups(group_file)
    else:
        groups = None
    return GlyTraitInputData(
        abundance_table=abundance_table, glycans=glycans, groups=groups
    )


@define(kw_only=True)
class GlyTraitInputData:
    """GlyTrait input data.

    Attributes:
        abundance_table: Abundance table as a pandas DataFrame.
        glycans: Glycans, either a dict of `Structure` objects or
            a dict of `Composition` objects.
        groups: Sample groups as a pandas Series.

    Notes:
        The glycan dict should have the same keys as the abundance table.
        The samples in the abundance table should have the same names as the
        samples in the groups.

    Raises:
        FileFormatError: If the abundance table has different samples as the groups,
            or if the glycan dict has different glycans as the abundance table.
    """

    abundance_table: AbundanceTable
    glycans: StructureDict | CompositionDict
    groups: Optional[GroupSeries] = None

    def __attrs_post_init__(self):
        if self.groups is not None:
            self._check_samples()
        self._check_glycans()

    def _check_samples(self):
        """Check if `groups` has the same samples as `abundance_table`."""
        if not set(self.abundance_table.index) == set(self.groups.index):
            msg = (
                "The samples in the abundance table and the groups should be the same."
            )
            raise FileFormatError(msg)

    def _check_glycans(self):
        """Check if `glycans` has the same glycans as `abundance_table`."""
        if not set(self.abundance_table.columns) == set(self.glycans.keys()):
            msg = (
                "The glycans in the abundance table and the glycans should be the same."
            )
            raise FileFormatError(msg)


# ===== Functions for loading individual files =====
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
    _check_columns("abundance table", filepath, ["Sample"])
    df = pd.read_csv(filepath, index_col="Sample")
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


# ===== Helper functions =====
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

    with open(filepath, encoding="utf-8-sig") as f:
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
    with open(filepath, encoding="utf-8-sig") as f:
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
