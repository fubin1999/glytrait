from collections.abc import Iterable

import numpy as np
import pandas as pd
from attrs import define, field
from numpy.typing import NDArray

from .glycan import NGlycan

valid_meta_properties = [
    ".",
    "isComplex",
    "isHighMannose",
    "isHybrid",
    "isBisecting",
    "is1Antennary",
    "is2Antennary",
    "is3Antennary",
    "is4Antennary",
    "totalAntenna",
    "coreFuc",
    "antennaryFuc",
    "totalFuc",
    "totalSia",
    "totalMan",
    "totalGal",
]


def build_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[NGlycan]
) -> pd.DataFrame:
    """Build a table of meta properties for glycans.

    The following meta properties are included:
        - isComplex: Whether the glycan is a complex type.
        - isHighMannose: Whether the glycan is a high-mannose type.
        - isHybrid: Whether the glycan is a hybrid type.
        - isBisecting: Whether the glycan has a bisection.
        - is1Antennary: Whether the glycan has a 1-antenna.
        - is2Antennary: Whether the glycan has a 2-antenna.
        - is3Antennary: Whether the glycan has a 3-antenna.
        - is4Antennary: Whether the glycan has a 4-antenna.
        - totalAntenna: The total number of antennae.
        - coreFuc: The number of fucoses on the core.
        - antennaryFuc: The number of fucoses on the antenna.
        - totalFuc: The total number of fucoses.
        - totalSia: The total number of sialic acids.
        - totalMan: The total number of mannoses.
        - totalGal: The total number of galactoses.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans. Could be the
            accession, or the composition.
        glycans (Iterable[NGlycan]): The glycans.

    Returns:
        pd.DataFrame: The table of meta properties, with glycan IDs as the index,
            and the property names as the columns.
    """
    meta_property_table = pd.DataFrame(index=list(glycan_ids))

    meta_property_table["isComplex"] = [g.is_complex() for g in glycans]
    meta_property_table["isHighMannose"] = [g.is_high_mannose() for g in glycans]
    meta_property_table["isHybrid"] = [g.is_hybrid() for g in glycans]
    meta_property_table["isBisecting"] = [g.is_bisecting() for g in glycans]
    meta_property_table["is1Antennary"] = [g.count_antenna() == 1 for g in glycans]
    meta_property_table["is2Antennary"] = [g.count_antenna() == 2 for g in glycans]
    meta_property_table["is3Antennary"] = [g.count_antenna() == 3 for g in glycans]
    meta_property_table["is4Antennary"] = [g.count_antenna() == 4 for g in glycans]
    meta_property_table["totalAntenna"] = [g.count_antenna() for g in glycans]
    meta_property_table["coreFuc"] = [g.count_core_fuc() for g in glycans]
    meta_property_table["antennaryFuc"] = [g.count_antennary_fuc() for g in glycans]
    meta_property_table["totalFuc"] = [g.count_fuc() for g in glycans]
    meta_property_table["totalSia"] = [g.count_sia() for g in glycans]
    meta_property_table["totalMan"] = [g.count_man() for g in glycans]
    meta_property_table["totalGal"] = [g.count_gal() for g in glycans]

    return meta_property_table


def _check_length(instance, attribute, value):
    if len(value) == 0:
        raise ValueError(f"`{attribute.name}` cannot be empty.")


def _check_meta_properties(instance, attribute, value):
    invalid_properties = set(value) - set(valid_meta_properties)
    if len(invalid_properties) > 0:
        raise ValueError(
            f"`{attribute.name}` contains invalid meta properties: "
            f"{', '.join(invalid_properties)}."
        )


_validators = [_check_length, _check_meta_properties]


@define
class TraitFormula:
    """The trait formula.

    Attributes:
        description (str): The description of the trait.
        name (str): The name of the trait.
        numerator_properties (tuple[str]): The meta properties in the numerator.
        denominator_properties (tuple[str]): The meta properties in the denominator.

    Examples:
        >>> formula = TraitFormula(
        ...     description="The ratio of high-mannose to complex glycans",
        ...     name="MHy",
        ...     numerator_properties=["isHighMannose"],
        ...     denominator_properties=["isComplex"],
        ... )
        >>> formula.initialize(meta_property_table)
        >>> trait_df[formula.name] = formula.calcu_trait(abundance_table)
    """

    description: str = field()
    name: str = field()
    numerator_properties: list[str] = field(converter=list, validator=_validators)
    denominator_properties: list[str] = field(converter=list, validator=_validators)
    _initialized = field(init=False, default=False)
    _numerator = field(init=False, default=None)
    _denominator = field(init=False, default=None)

    def initialize(self, meta_property_table: pd.DataFrame) -> None:
        """Initialize the trait formula.

        Args:
            meta_property_table (pd.DataFrame): The table of meta properties generated
                by `build_meta_property_table`.
        """
        self._numerator = self._initialize(
            meta_property_table, self.numerator_properties
        )
        self._denominator = self._initialize(
            meta_property_table, self.denominator_properties
        )
        self._initialized = True

    @staticmethod
    def _initialize(
        meta_property_table: pd.DataFrame, properties: list[str]
    ) -> NDArray:
        if len(properties) == 1 and properties[0] == ".":
            return np.ones_like(meta_property_table.index)
        else:
            return meta_property_table[properties].prod(axis=1)

    def calcu_trait(self, abundance_table: pd.DataFrame) -> NDArray:
        """Calculate the trait.

        Args:
            abundance_table (pd.DataFrame): The glycan abundance table, with samples as index,
                and glycans as columns.

        Returns:
            NDArray: An array of trait values for each sample.
        """
        if not self._initialized:
            raise RuntimeError("TraitFormula is not initialized.")

        numerator = abundance_table.values @ self._numerator
        denominator = abundance_table.values @ self._denominator
        return numerator / denominator
