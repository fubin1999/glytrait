"""Functions for loading data for GlyTrait.

This module provides functions for loading data for GlyTrait.

Classes:
    GlyTraitInputData: Encapsulates all the input data for GlyTrait.

Functions:
    load_input_data: Load all the input data for GlyTrait, including
        abundance table, glycans, and groups.
        Returns a `GlyTraitInputData` object.
"""

from __future__ import annotations

from collections.abc import Iterable, Callable
from typing import Optional, Literal, Protocol, cast, ClassVar

import pandas as pd
from attrs import define, field
from numpy import dtype

from glytrait.data_type import (
    AbundanceTable,
    GroupSeries,
)
from glytrait.exception import DataInputError
from glytrait.glycan import (
    parse_structures,
    parse_compositions,
    StructureDict,
    CompositionDict,
    GlycanDict,
)

__all__ = [
    "GlyTraitInputData",
    "AbundanceCSVLoader",
    "GroupsCSVLoader",
    "GlycanCSVLoader",
    "load_input_data",
]


# ===== Loaders =====
class AbundanceLoaderProto(Protocol):
    def load(self) -> AbundanceTable:
        ...


class GlycanLoaderProto(Protocol):
    def load(self) -> StructureDict | CompositionDict:
        ...


class GroupsLoaderProto(Protocol):
    def load(self) -> GroupSeries:
        ...


@define
class DFValidator:
    """Validator for pandas DataFrame.

    Attributes:
        must_have: List of column names that must be in the DataFrame.
        unique: List of column names that must be unique in the DataFrame.
        types: Dictionary of column names and their expected types.
        default_type: Default type for columns not in `types`.
    """

    must_have: list[str] = field(kw_only=True, factory=list)
    unique: list[str] = field(kw_only=True, factory=list)
    types: dict[str, str] = field(kw_only=True, factory=dict)
    default_type: dtype | str = field(kw_only=True, default=None)

    def __call__(self, df: pd.DataFrame) -> None:
        """Validate the DataFrame.

        Raises:
            DataInputError: If one of the following conditions is met:
            - The DataFrame does not have all the columns specified in `must_have`.
        """
        self._test_must_have_columns(df)
        self._test_unique_columns(df)
        self._test_type_check(df)
        self._test_default_type(df)

    def _test_must_have_columns(self, df: pd.DataFrame):
        if missing := {col for col in self.must_have if col not in df.columns}:
            msg = f"The following columns are missing: {', '.join(missing)}."
            raise DataInputError(msg)

    def _test_unique_columns(self, df: pd.DataFrame):
        if non_unique := [
            col for col in self.unique if col in df and df[col].duplicated().any()
        ]:
            msg = f"The following columns are not unique: {', '.join(non_unique)}."
            raise DataInputError(msg)

    def _test_type_check(self, df: pd.DataFrame):
        if wrong_type_cols := [
            col
            for col, dtype in self.types.items()
            if col in df and df[col].dtype != dtype
        ]:
            expected = {col: self.types[col] for col in wrong_type_cols}
            got = {col: df[col].dtype for col in wrong_type_cols}
            msg = (
                f"The following columns have incorrect types: {', '.join(wrong_type_cols)}. "
                f"Expected types: {expected}, got: {got}."
            )
            raise DataInputError(msg)

    def _test_default_type(self, df: pd.DataFrame):
        if self.default_type is None:
            return
        cols_to_check = set(df.columns) - set(self.types)
        if wrong_type_cols := [
            col for col in cols_to_check if df[col].dtype != self.default_type
        ]:
            got = {col: df[col].dtype for col in wrong_type_cols}
            msg = (
                f"The following columns have incorrect types: {', '.join(wrong_type_cols)}. "
                f"Expected types: {self.default_type}, got: {got}."
            )
            raise DataInputError(msg)


@define
class AbundanceCSVLoader:
    """Loader for abundance table from a csv file."""

    filepath: str
    validator: DFValidator = field(
        kw_only=True,
        default=DFValidator(
            must_have=["Sample"],
            unique=["Sample"],
            types={"Sample": "object"},
            default_type=dtype("float64"),
        ),
    )

    def load(self) -> AbundanceTable:
        """Returns the abundance table loaded from the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            FileFormatError: If the file format is incorrect.
                This includes: (1) the file is empty; (2) the file could not be parsed;
                (3) the "Sample" column is not found.
        """
        df = _read_csv(self.filepath)
        self.validator(df)
        return AbundanceTable(df.set_index("Sample"))


@define
class GroupsCSVLoader:
    """Loader for groups from a csv file."""

    filepath: str
    validator: DFValidator = field(
        kw_only=True,
        default=DFValidator(
            must_have=["Group", "Sample"],
            unique=["Sample"],
            types={"Sample": "object"},
        ),
    )

    def load(self) -> GroupSeries:
        """Returns the groups loaded from the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            FileFormatError: If the file format is incorrect.
                This includes: (1) the file is empty; (2) the file could not be parsed;
                (3) the "Sample" column is not found; (4) the "Group" column is not found.
        """
        df = _read_csv(self.filepath)
        self.validator(df)
        return GroupSeries(df.set_index("Sample")["Group"])


GlycanParserType = Callable[
    [Iterable[tuple[str, str]]], StructureDict | CompositionDict
]


@define
class GlycanCSVLoader:
    """Loader for structures or compositions from a csv file."""

    filepath: str
    mode: Literal["structure", "composition"] = field(kw_only=True)
    validator: DFValidator = field(kw_only=True, default=None)
    parser: GlycanParserType = field(kw_only=True, default=None)

    def __attrs_post_init__(self):
        if self.parser is None:
            self.parser = self._glycan_parser_factory(self.mode)
        if self.validator is None:
            self.validator = self._validator_factory(self.mode)

    @staticmethod
    def _glycan_parser_factory(
        mode: Literal["structure", "composition"]
    ) -> GlycanParserType:
        parser = parse_structures if mode == "structure" else parse_compositions
        return cast(GlycanParserType, parser)

    @staticmethod
    def _validator_factory(mode: Literal["structure", "composition"]) -> DFValidator:
        glycan_col = "Structure" if mode == "structure" else "Composition"
        return DFValidator(
            must_have=["GlycanID", glycan_col],
            unique=["GlycanID", glycan_col],
            types={"GlycanID": "object", glycan_col: "object"},
        )

    def load(self) -> GlycanDict:
        """Returns the glycans loaded from the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            FileFormatError: If the file format is incorrect.
                This includes: (1) the file is empty; (2) the file could not be parsed;
                (3) the "GlycanID" column is not found; (4) the "Structure" or
                "Composition" column is not found.
            StructureParseError: If any structure cannot be parsed,
                when `mode` is "structure".
            CompositionParseError: If any composition cannot be parsed,
                when `mode` is "composition".
        """
        df = _read_csv(self.filepath)
        self.validator(df)
        ids = df["GlycanID"].to_list()
        glycan_col = self.mode.capitalize()
        try:
            strings = df[glycan_col].to_list()
        except KeyError as e:
            raise DataInputError(f"The '{glycan_col}' column is not found.") from e
        return self.parser(zip(ids, strings))


def _read_csv(filepath: str) -> pd.DataFrame:
    try:
        return pd.read_csv(filepath)
    # FileNotFoundError is not caught here
    except pd.errors.EmptyDataError as e:
        raise DataInputError("Empty CSV file.") from e
    except pd.errors.ParserError as e:
        raise DataInputError("This CSV file could not be parsed.") from e


# ===== Input data =====
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


# ===== Input data validators =====
ValidatorType = Callable[[GlyTraitInputData], None]


@define
class InputDataValidator:
    """Validator for the input data.

    This class only validates the interaction between the input data.
    The format of each input data is validated by the loaders.
    """

    validators: ClassVar[list[ValidatorType]] = []

    def __call__(self, input_data: GlyTraitInputData) -> None:
        for validator in self.validators:
            validator(input_data)

    @classmethod
    def add_validator(cls, validator: ValidatorType) -> ValidatorType:
        """A decorator to add a validator to the list of validators."""
        cls.validators.append(validator)
        return validator


@InputDataValidator.add_validator
def same_samples_in_abundance_and_groups(input_data: GlyTraitInputData) -> None:
    """Check if the abundance table and the groups have the same samples.

    Raises:
        DataInputError: If the samples in the abundance table and the groups are different.
    """
    if input_data.groups is None:
        return
    abund_samples = set(input_data.abundance_table.index)
    groups_samples = set(input_data.groups.index)
    if abund_samples != groups_samples:
        samples_in_abund_not_in_groups = abund_samples - groups_samples
        samples_in_groups_not_in_abund = groups_samples - abund_samples
        msg = ""
        if samples_in_abund_not_in_groups:
            msg += (
                f"The following samples are in the abundance table but not in the groups: "
                f"{', '.join(samples_in_abund_not_in_groups)}. "
            )
        if samples_in_groups_not_in_abund:
            msg += (
                f"The following samples are in the groups but not in the abundance table: "
                f"{', '.join(samples_in_groups_not_in_abund)}."
            )
        raise DataInputError(msg)


@InputDataValidator.add_validator
def all_glycans_have_structures_or_compositions(input_data: GlyTraitInputData) -> None:
    """Check if all glycans in the abundance table have structures or compositions.

    Glycans in the structure or composition dict that are not in the abundance table
    are not checked.

    Raises:
        DataInputError: If any glycan in the abundance table does not
            have a structure or composition.
    """
    abund_glycans = set(input_data.abundance_table.columns)
    glycans = set(input_data.glycans.keys())
    if diff := abund_glycans - glycans:
        msg = (
            f"The following glycans in the abundance table do not have structures or "
            f"compositions: {', '.join(diff)}."
        )
        raise DataInputError(msg)


# ===== The highest-level API =====
def load_input_data(
    *,
    abundance_file: str,
    glycan_file: str,
    group_file: Optional[str] = None,
    abundance_loader: Optional[AbundanceLoaderProto] = None,
    glycan_loader: Optional[GlycanLoaderProto] = None,
    group_loader: Optional[GroupsLoaderProto] = None,
    mode: Literal["structure", "composition"] = "structure",
    validator: ValidatorType = InputDataValidator(),
) -> GlyTraitInputData:
    """Load all the input data for GlyTrait.

    Notes:
        - If `group_file` is not provided, the `groups` attribute of the returned
        `GlyTraitInputData` will be `None`.
        - If `abundance_loader`, `glycan_loader`, or `group_loader` are not provided,
        the default CSV loaders will be used.
        - If `abundance_loader`, `glycan_loader`, or `group_loader` are provided,
        the three "file" arguments will be ignored.
        - The three "file" arguments will be removed in the future.
        It exists for backward compatibility.
        To use the new API, provide the loaders directly, and pass "" as the file paths.

    Args:
        abundance_file: Path to the abundance table file.
        glycan_file: Path to the glycans file.
        group_file: Path to the groups file. Optional.
        abundance_loader: Loader for the abundance table.
        glycan_loader: Loader for the glycans.
        group_loader: Loader for the groups. Optional.
        mode: Either "structure" or "composition".
        validator: Validator for the input data.

    Returns:
        GlyTraitInputData: Input data for GlyTrait.
    """
    if abundance_loader is None:
        abundance_loader = AbundanceCSVLoader(filepath=abundance_file)
    if glycan_loader is None:
        glycan_loader = GlycanCSVLoader(filepath=glycan_file, mode=mode)
    if group_loader is None and group_file is not None:
        group_loader = GroupsCSVLoader(filepath=group_file)

    abundance_table = abundance_loader.load()
    glycans = glycan_loader.load()
    groups = group_loader.load() if group_loader else None

    input_data = GlyTraitInputData(
        abundance_table=abundance_table,
        glycans=glycans,
        groups=groups,
    )

    validator(input_data)
    return input_data
