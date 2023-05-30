import itertools
import re
from collections.abc import Iterable, Generator
from importlib.resources import files
from pathlib import Path
from typing import Optional

import attrs.validators
import numpy as np
import pandas as pd
from attrs import define, field
from numpy.typing import NDArray

from .exception import FormulaError
from .glycan import NGlycan

DEFAULT_FORMULA_FILE = "trait_forluma.txt"
FORMULA_TEMPLATE_FILE = "trait_formula_template.txt"

basic_meta_properties = [
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
    "hasAntennaryFuc",
    "totalFuc",
    "hasFuc",
    "noFuc",
    "totalSia",
    "hasSia",
    "noSia",
    "totalMan",
    "totalGal",
]

sia_linkage_meta_properties = [
    "a23Sia",
    "a26Sia",
    "hasa23Sia",
    "hasa26Sia",
    "noa23Sia",
    "noa26Sia",
]

meta_properties = basic_meta_properties + sia_linkage_meta_properties


def build_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[NGlycan], sia_linkage: bool = False
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
        - hasAntennaryFuc: Whether the glycan has any fucoses on the antenna.
        - totalFuc: The total number of fucoses.
        - hasFuc: Whether the glycan has any fucoses.
        - noFuc: Whether the glycan has no fucoses.
        - totalSia: The total number of sialic acids.
        - hasSia: Whether the glycan has any sialic acids.
        - noSia: Whether the glycan has no sialic acids.
        - totalMan: The total number of mannoses.
        - totalGal: The total number of galactoses.

    If `sia_linkage` is True, the following meta properties are also included:
        - a23Sia: The number of sialic acids with an alpha-2,3 linkage.
        - a26Sia: The number of sialic acids with an alpha-2,6 linkage.
        - hasa23Sia: Whether the glycan has any sialic acids with an alpha-2,3 linkage.
        - hasa26Sia: Whether the glycan has any sialic acids with an alpha-2,6 linkage.
        - noa23Sia: Whether the glycan has no sialic acids with an alpha-2,3 linkage.
        - noa26Sia: Whether the glycan has no sialic acids with an alpha-2,6 linkage.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans. Could be the
            accession, or the composition.
        glycans (Iterable[NGlycan]): The glycans.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

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
    meta_property_table["hasAntennaryFuc"] = [
        g.count_antennary_fuc() > 0 for g in glycans
    ]
    meta_property_table["totalFuc"] = [g.count_fuc() for g in glycans]
    meta_property_table["hasFuc"] = [g.count_fuc() > 0 for g in glycans]
    meta_property_table["noFuc"] = [g.count_fuc() == 0 for g in glycans]
    meta_property_table["totalSia"] = [g.count_sia() for g in glycans]
    meta_property_table["hasSia"] = [g.count_sia() > 0 for g in glycans]
    meta_property_table["noSia"] = [g.count_sia() == 0 for g in glycans]
    meta_property_table["totalMan"] = [g.count_man() for g in glycans]
    meta_property_table["totalGal"] = [g.count_gal() for g in glycans]

    if sia_linkage:
        meta_property_table["a23Sia"] = [g.count_a23_sia() for g in glycans]
        meta_property_table["a26Sia"] = [g.count_a26_sia() for g in glycans]
        meta_property_table["hasa23Sia"] = [g.count_a23_sia() > 0 for g in glycans]
        meta_property_table["hasa26Sia"] = [g.count_a26_sia() > 0 for g in glycans]
        meta_property_table["noa23Sia"] = [g.count_a23_sia() == 0 for g in glycans]
        meta_property_table["noa26Sia"] = [g.count_a26_sia() == 0 for g in glycans]

    return meta_property_table


def _check_length(instance, attribute, value):
    if len(value) == 0:
        raise FormulaError(f"`{attribute.name}` cannot be empty.")


def _check_meta_properties(instance, attribute, value):
    invalid_properties = set(value) - set(meta_properties)
    if len(invalid_properties) > 0:
        raise FormulaError(
            f"`{attribute.name}` contains invalid meta properties: "
            f"{', '.join(invalid_properties)}."
        )


def _check_numerator(instance, attribute, value):
    if "." in value:
        raise FormulaError("'.' should not be used in the numerator.")


def _check_denominator(instance, attribute, value):
    if "." in value and len(value) > 1:
        raise FormulaError(
            "'.' should not be used with other meta properties in the denominator."
        )


@define
class TraitFormula:
    """The trait formula.

    Attributes:
        description (str): The description of the trait.
        name (str): The name of the trait.
        numerator_properties (tuple[str]): The meta properties in the numerator.
        denominator_properties (tuple[str]): The meta properties in the denominator.
        coefficient (float): The coefficient of the trait.

    Examples:
        >>> from glytrait.trait import TraitFormula
        >>> formula = TraitFormula(
        ...     description="The ratio of high-mannose to complex glycans",
        ...     name="MHy",
        ...     numerator_properties=["isHighMannose"],
        ...     denominator_properties=["isComplex"],
        ... )
        >>> formula.initialize(meta_property_table)
        >>> trait_df[formula.name] = glytrait.trait.calcu_trait(abundance_table)
    """

    description: str = field()
    name: str = field()
    numerator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_numerator],
    )
    denominator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_denominator],
    )
    coefficient: float = field(default=1.0, validator=attrs.validators.gt(0))
    _sia_linkage: bool = field(init=False, default=False)
    _initialized = field(init=False, default=False)
    _numerator = field(init=False, default=None)
    _denominator = field(init=False, default=None)

    def __attrs_post_init__(self):
        for prop in itertools.chain(
            self.numerator_properties, self.denominator_properties
        ):
            if prop in sia_linkage_meta_properties:
                self._sia_linkage = True
                break

    @property
    def sia_linkage(self) -> bool:
        """Whether the formula contains sia linkage meta properties."""
        return self._sia_linkage

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
        if np.any(denominator == 0):
            return np.zeros_like(numerator)
        return numerator / denominator * self.coefficient


def _load_formulas(
    formula_file_reader: Iterable[str],
) -> Generator[TraitFormula, None, None]:
    """Load the formulas from a file.

    Args:
        formula_file_reader (Iterable[str]): The path of the formula file.

    Returns:
        Generator[TraitFormula, None, None]: The generator of the formulas.
    """
    description = None
    expression = None
    for line in formula_file_reader:
        if line.startswith("@"):
            if description is not None:
                raise FormulaError(
                    f"Invalid line: '{line}'.\n"
                    f"One description line must follow a formula line."
                )
            description = line[1:].strip()

        if line.startswith("$"):
            if description is None:
                raise FormulaError(
                    f"Invalid line: '{line}'.\n"
                    "One formula line must follow a description line."
                )
            expression = line[1:].strip()
            name, num_prop, den_prop, coef = _parse_expression(expression)
            try:
                yield TraitFormula(
                    description=description,
                    name=name,
                    numerator_properties=num_prop,
                    denominator_properties=den_prop,
                    coefficient=coef,
                )
            except FormulaError as e:
                raise FormulaError(
                    f"Invalid line: '{line}'.\n"
                    f"Error in formula '{expression}': {e}"
                )
            description = None
            expression = None

    if description is not None and expression is None:
        raise FormulaError(
            f"Invalid line: '{line}'.\n"
            f"One description line must follow a formula line."
        )


def load_default_formulas() -> Generator[TraitFormula, None, None]:
    """Load the default formulas."""
    file_reader = files("glytrait").joinpath(DEFAULT_FORMULA_FILE).open("r")
    yield from _load_formulas(file_reader)


def load_user_formulas(file: str) -> Generator[TraitFormula, None, None]:
    """Load the user-defined formulas.

    If the formula name is duplicated, the first one will be used.

    Args:
        file (str): The path of the user-defined formula file.

    Yields:
        The formulas parsed.
    """
    formulas_parsed: set[str] = set()
    file_reader = Path(file).open("r")
    for formula in _load_formulas(file_reader):
        if formula.name in formulas_parsed:
            continue
        yield formula
        formulas_parsed.add(formula.name)


def load_formulas(user_file: Optional[str] = None) -> list[TraitFormula]:
    """Load both the default formulas and the user-defined formulas."""
    default_formulas = list(load_default_formulas())
    if user_file is None:
        return default_formulas

    default_formula_names = [f.name for f in default_formulas]
    user_formulas = list(load_user_formulas(user_file))

    formulas: list[TraitFormula] = []
    formulas.extend(default_formulas)
    for f in user_formulas:
        if f.name in default_formula_names:
            continue
        formulas.append(f)
    return formulas


def save_trait_formula_template(dirpath: str) -> None:
    """Copy the template of trait formula file to the given path.

    Args:
        dirpath (str): The path to save the template.
    """
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    file = dirpath / "trait_formula.txt"
    content = files("glytrait").joinpath(FORMULA_TEMPLATE_FILE).open("r").read()
    with open(file, "w") as f:
        f.write(content)


def _parse_expression(expr: str) -> tuple[str, list[str], list[str], float]:
    """Parse the expression of a formula.

    Args:
        expr (str): The expression of a formula.

    Returns:
        tuple[str, list[str], list[str], float]: The name, numerator properties, denominator
            properties, and the coefficient of the formula.
    """
    if "//" in expr:
        pattern = r"(\w+) = \((.+)\) // \((.+)\)"  # Expression with the "//" shortcut
    else:
        pattern = r"(\w+) = \((.+)\) / \((.+)\)"  # Normal expression

    match = re.match(pattern, expr)
    if match is None:
        raise FormulaError(f"Invalid expression: '{expr}'")
    name, num_prop, den_prop = match.groups()

    num_prop = num_prop.split("*")
    num_prop = [p.strip() for p in num_prop]
    den_prop = den_prop.split("*")
    den_prop = [p.strip() for p in den_prop]

    # Check if there are invalid characters in the properties
    for prop in num_prop + den_prop:
        if re.search(r"\W", prop) and prop != ".":
            raise FormulaError(f"Invalid expression: '{expr}'")

    # If "//" is used, we need to add the denominator properties to the numerator properties.
    if "//" in expr:
        num_prop.extend(den_prop)

    # Parse the coefficient
    if ") *" in expr:
        coef_pattern = r"\) \* (\d+/\d+|\d+(\.\d+)?)"
        match = re.search(coef_pattern, expr)
        if match is None:
            raise FormulaError(f"Invalid expression: '{expr}'")
        coef = eval(match.group(1))
    else:
        coef = 1.0

    return name, num_prop, den_prop, coef


def calcu_trait(
    abund_df: pd.DataFrame,
    meta_prop_df: pd.DataFrame,
    formulas: list[TraitFormula],
) -> pd.DataFrame:
    """Calculate the trait values.

    Args:
        abund_df (pd.DataFrame): The abundance table, with samples as index and glycan IDs
            as columns.
        meta_prop_df (pd.DataFrame): The table of meta properties generated by
            `build_meta_property_table`.
        formulas (list[TraitFormula]): The trait formulas.

    Returns:
        pd.DataFrame: The trait values, with samples as index and trait names as columns.
    """
    trait_series: list[pd.Series] = []
    for formula in formulas:
        formula.initialize(meta_prop_df)
        trait_s = pd.Series(
            data=formula.calcu_trait(abund_df), index=abund_df.index, name=formula.name
        )
        trait_series.append(trait_s)
    return pd.concat(trait_series, axis=1)
