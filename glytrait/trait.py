import itertools
import re
from collections.abc import Iterable, Generator
from importlib.resources import files
from pathlib import Path
from typing import Optional, Literal

import attrs.validators
import numpy as np
import pandas as pd
from attrs import define, field
from numpy.typing import NDArray

from .exception import FormulaError
from .glycan import NGlycan, Composition

default_struc_formula_file = files("glytrait.resources").joinpath("struc_formula.txt")
default_comp_formula_file = files("glytrait.resources").joinpath("comp_formula.txt")
formula_template_file = files("glytrait.resources").joinpath(
    "trait_formula_template.txt"
)

basic_struc_meta_properties = {
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
}

sia_struc_meta_properties = {
    "a23Sia",
    "a26Sia",
    "hasa23Sia",
    "hasa26Sia",
    "noa23Sia",
    "noa26Sia",
}

struc_meta_properties = basic_struc_meta_properties | sia_struc_meta_properties

basic_comp_meta_properties = {
    ".",
    "isHighBranching",
    "isLowBranching",
    "totalSia",
    "totalFuc",
    "totalGal",
    "hasSia",
    "hasFuc",
    "hasGal",
    "noSia",
    "noFuc",
    "noGal",
}

sia_comp_meta_properties = {
    "a23Sia",
    "a26Sia",
    "hasa23Sia",
    "hasa26Sia",
    "noa23Sia",
    "noa26Sia",
}

comp_meta_properties = basic_comp_meta_properties | sia_comp_meta_properties


def build_meta_property_table(
    glycan_ids: Iterable[str],
    glycans: Iterable[NGlycan | Composition],
    mode: Literal["composition", "structure"],
    sia_linkage: bool = False,
) -> pd.DataFrame:
    """Build a table of meta properties for glycans.

    Args:
        glycan_ids (Iterable[str]): The IDs of the glycans.
        glycans (Iterable[NGlycan]): The glycans.
        mode (Literal["composition", "structure"]): The calculation mode.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

    Returns:
        pd.DataFrame: The table of meta properties, with Compositions as the index,
    """
    if mode == "composition":
        return _build_comp_meta_property_table(glycan_ids, glycans, sia_linkage)
    elif mode == "structure":
        return _build_struc_meta_property_table(glycan_ids, glycans, sia_linkage)
    else:
        raise ValueError(f"Invalid type: {mode}")


def _build_struc_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[NGlycan], sia_linkage: bool = False
) -> pd.DataFrame:
    """Build a table of structural meta properties for glycans.

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
        pd.DataFrame: The table of meta properties, with Compositions as the index,
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


def _build_comp_meta_property_table(
    glycan_ids: Iterable[str], glycans: Iterable[Composition], sia_linkage: bool = False
) -> pd.DataFrame:
    """Build a table of compositional meta properties for glycans.

    The following meta properties are included:
        - isHighBranching: Whether the glycan has a high branching.
        - isLowBranching: Whether the glycan has a low branching.
        - totalSia: The total number of sialic acids.
        - totalFuc: The total number of fucoses.
        - totalGal: The total number of galactoses.
        - hasSia: Whether the glycan has any sialic acids.
        - hasFuc: Whether the glycan has any fucoses.
        - hasGal: Whether the glycan has any galactoses.
        - noSia: Whether the glycan has no sialic acids.
        - noFuc: Whether the glycan has no fucoses.
        - noGal: Whether the glycan has no galactoses.

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
        glycans (Iterable[Composition]): The glycans.
        sia_linkage (bool, optional): Whether to include the sialic acid linkage
            meta properties. Defaults to False.

    Returns:
        pd.DataFrame: The table of meta properties, with Compositions as the index,
            and the property names as the columns.
    """
    meta_property_table = pd.DataFrame(index=list(glycan_ids))

    meta_property_table["isHighBranching"] = [g.is_high_branching() for g in glycans]
    meta_property_table["isLowBranching"] = [g.is_low_branching() for g in glycans]
    meta_property_table["totalSia"] = [g.count_sia() for g in glycans]
    meta_property_table["totalFuc"] = [g.count_fuc() for g in glycans]
    meta_property_table["totalGal"] = [g.count_gal() for g in glycans]
    meta_property_table["hasSia"] = [g.count_sia() > 0 for g in glycans]
    meta_property_table["hasFuc"] = [g.count_fuc() > 0 for g in glycans]
    meta_property_table["hasGal"] = [g.count_gal() > 0 for g in glycans]
    meta_property_table["noSia"] = [g.count_sia() == 0 for g in glycans]
    meta_property_table["noFuc"] = [g.count_fuc() == 0 for g in glycans]
    meta_property_table["noGal"] = [g.count_gal() == 0 for g in glycans]

    if sia_linkage:
        meta_property_table["a23Sia"] = [g.count_a23_sia() for g in glycans]
        meta_property_table["a26Sia"] = [g.count_a26_sia() for g in glycans]
        meta_property_table["hasa23Sia"] = [g.count_a23_sia() > 0 for g in glycans]
        meta_property_table["hasa26Sia"] = [g.count_a26_sia() > 0 for g in glycans]
        meta_property_table["noa23Sia"] = [g.count_a23_sia() == 0 for g in glycans]
        meta_property_table["noa26Sia"] = [g.count_a26_sia() == 0 for g in glycans]

    return meta_property_table


def _check_length(instance, attribute, value):
    """Validator for `TraitFormula`."""
    if len(value) == 0:
        raise FormulaError(f"`{attribute.name}` cannot be empty.")


def _check_meta_properties(instance, attribute, value):
    """Validator for `TraitFormula`."""
    if instance.type == "composition":
        valid_meta_properties = comp_meta_properties
    else:
        valid_meta_properties = struc_meta_properties
    invalid_properties = set(value) - set(valid_meta_properties)
    if len(invalid_properties) > 0:
        raise FormulaError(
            f"`{attribute.name}` contains invalid meta properties: "
            f"{', '.join(invalid_properties)}."
        )


def _check_numerator(instance, attribute, value):
    """Validator for `TraitFormula`."""
    if "." in value:
        raise FormulaError("'.' should not be used in the numerator.")


def _check_denominator(instance, attribute, value):
    """Validator for `TraitFormula`."""
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
        type (str): The type of the trait. Either "structure" or "composition".
        numerator_properties (tuple[str]): The meta properties in the numerator.
        denominator_properties (tuple[str]): The meta properties in the denominator.
        coefficient (float): The coefficient of the trait.

    Examples:
        >>> from glytrait.trait import TraitFormula
        >>> formula = TraitFormula(
        ...     description="The ratio of high-mannose to complex glycans",
        ...     name="MHy",
        ...     type="structure",
        ...     numerator_properties=["isHighMannose"],
        ...     denominator_properties=["isComplex"],
        ... )
        >>> formula.initialize(meta_property_table)
        >>> trait_df[formula.name] = glytrait.trait.calcu_derived_trait(abundance_table)
    """

    description: str = field()
    name: str = field()
    type: str = field(validator=attrs.validators.in_(["structure", "composition"]))
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
        self._sia_linkage = self._init_sia_linkage()

    def _init_sia_linkage(self) -> bool:
        """Whether the formula contains sia linkage meta properties."""
        if self.type == "composition":
            sia_meta_properties = sia_comp_meta_properties
        else:
            sia_meta_properties = sia_struc_meta_properties
        for prop in itertools.chain(
            self.numerator_properties, self.denominator_properties
        ):
            if prop in sia_meta_properties:
                return True
        return False

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
            return np.array([np.nan] * len(abundance_table))
        return numerator / denominator * self.coefficient


def load_formulas(
    type: Literal["structure", "composition"], user_file: Optional[str] = None
) -> list[TraitFormula]:
    """Load both the default formulas and the user-defined formulas.

    Args:
        type (Literal["structure", "composition"]): The type of the formulas.
        user_file (Optional[str], optional): The path of the user-defined formula file.

    Returns:
        list[TraitFormula]: The formulas.
    """
    default_formulas = list(_load_default_formulas(type=type))
    if user_file is None:
        return default_formulas

    default_formula_names = [f.name for f in default_formulas]
    user_formulas = list(_load_user_formulas(user_file, type=type))

    formulas: list[TraitFormula] = []
    formulas.extend(default_formulas)
    for f in user_formulas:
        if f.name in default_formula_names:
            continue
        formulas.append(f)
    return formulas


def _load_default_formulas(
    type: Literal["structure", "composition"]
) -> Generator[TraitFormula, None, None]:
    """Load the default formulas.

    Args:
        type (Literal["structure", "composition"]): The type of the formulas.

    Yields:
        The formulas parsed.
    """
    if type == "composition":
        default_file = default_comp_formula_file
    else:
        default_file = default_struc_formula_file
    file_reader = default_file.open("r")
    yield from _load_formulas(file_reader, type=type)


def _load_user_formulas(
    file: str, type: Literal["structure", "composition"]
) -> Generator[TraitFormula, None, None]:
    """Load the user-defined formulas.

    If the formula name is duplicated, the first one will be used.

    Args:
        file (str): The path of the user-defined formula file.
        type (Literal["structure", "composition"]): The type of the formulas.

    Yields:
        The formulas parsed.
    """
    formulas_parsed: set[str] = set()
    file_reader = Path(file).open("r")
    for formula in _load_formulas(file_reader, type=type):
        if formula.name in formulas_parsed:
            continue
        yield formula
        formulas_parsed.add(formula.name)


def _load_formulas(
    formula_file_reader: Iterable[str], type: Literal["structure", "composition"]
) -> Generator[TraitFormula, None, None]:
    """Load the formulas from a file.

    Args:
        formula_file_reader (Iterable[str]): The path of the formula file.
        type (Literal["structure", "composition"]): The type of the formulas.

    Returns:
        Generator[TraitFormula, None, None]: The generator of the formulas.
    """
    description = None
    expression = None
    for line in formula_file_reader:
        line = line.strip()
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
                    type=type,
                    numerator_properties=num_prop,
                    denominator_properties=den_prop,
                    coefficient=coef,
                )
            except FormulaError as e:
                raise FormulaError(
                    f"Invalid line: '{line}'.\n" f"Error in formula '{expression}': {e}"
                )
            description = None
            expression = None

    if description is not None and expression is None:
        raise FormulaError(
            f"Invalid line: '{line}'.\n"
            f"One description line must follow a formula line."
        )


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


def save_trait_formula_template(dirpath: str) -> None:
    """Copy the template of trait formula file to the given path.

    Args:
        dirpath (str): The path to save the template.
    """
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    file = dirpath / "trait_formula.txt"
    content = formula_template_file.open("r").read()
    with open(file, "w", encoding="utf8") as f:
        f.write(content)


def save_builtin_formula(dirpath: str) -> None:
    """Copy the builtin formula file to the given path.

    Args:
        dirpath (str): The path to save the built-in formula file.
    """
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    struc_file = dirpath / "struc_builtin_formulas.txt"
    comp_file = dirpath / "comp_builtin_formulas.txt"
    struc_content = default_struc_formula_file.open("r").read()
    comp_content = default_comp_formula_file.open("r").read()
    with open(struc_file, "w", encoding="utf8") as f:
        f.write(struc_content)
    with open(comp_file, "w", encoding="utf8") as f:
        f.write(comp_content)


def calcu_derived_trait(
    abund_df: pd.DataFrame,
    meta_prop_df: pd.DataFrame,
    formulas: list[TraitFormula],
) -> pd.DataFrame:
    """Calculate the derived trait values.

    Args:
        abund_df (pd.DataFrame): The abundance table, with samples as index and Compositions
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
            data=formula.calcu_trait(abund_df),
            index=abund_df.index,
            name=formula.name,
            dtype=float,
        )
        trait_series.append(trait_s)
    derived_trait_df = pd.concat(trait_series, axis=1)
    derived_trait_df = derived_trait_df.round(6)
    return derived_trait_df


def filter_derived_trait(trait_df: pd.DataFrame) -> pd.DataFrame:
    """Rule out the invalid traits.

    A trait is invalid if it:
    1. Has the same value for all samples.
    2. Is NaN for all samples.

    Args:
        trait_df (pd.DataFrame): The trait values, with samples as index and trait names
            as columns.

    Returns:
        pd.DataFrame: The filtered trait values.
    """
    trait_df = _filter_all_same(trait_df)
    trait_df = _filter_all_nan(trait_df)
    return trait_df


def _filter_all_same(trait_df: pd.DataFrame) -> pd.DataFrame:
    """Rule out the traits that have the same value for all samples."""
    return trait_df.loc[:, trait_df.nunique() != 1]


def _filter_all_nan(trait_df: pd.DataFrame) -> pd.DataFrame:
    """Rule out the traits that are NaN for all samples."""
    return trait_df.loc[:, trait_df.notna().any()]
