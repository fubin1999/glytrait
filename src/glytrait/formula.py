"""This module defines the `TraitFormula` class, and other related functions.

The `TraitFormula` class is used to represent a trait formula.
It is the core of calculating derived traits.

Classes:
    TraitFormula: The trait formula.

Functions:
    load_formulas: Load all trait formulas from default formula files and custom formulas files.
    save_trait_formulas_tepmlate: Save a template of trait formulas to a file.
    save_builtin_formula: Save a built-in trait formulas to a file.
"""
from __future__ import annotations

import itertools
import re
from importlib.resources import files
from pathlib import Path
from typing import Literal, Optional, Generator, Iterable, Annotated

import attrs
import numpy as np
import pandas as pd
from attrs import field, frozen
from numpy.typing import NDArray

import glytrait
from glytrait.exception import FormulaError
from glytrait.meta_property import available_meta_properties

default_struc_formula_file = files("glytrait.resources").joinpath("struc_formula.txt")
default_comp_formula_file = files("glytrait.resources").joinpath("comp_formula.txt")
formula_template_file = files("glytrait.resources").joinpath(
    "trait_formula_template.txt"
)


def _check_length(instance, attribute, value):
    """Validator for `TraitFormula`."""
    if len(value) == 0:
        # The attribute name starts with "_", so we need to remove it.
        raise FormulaError(f"`{attribute.name[1:]}` cannot be empty.")


def _check_meta_properties(instance: TraitFormula, attribute, value):
    """Validator for `TraitFormula`."""
    valid_meta_properties = available_meta_properties(
        instance.type, sia_linkage=True  # type: ignore
    )
    invalid_properties = set(value) - set(valid_meta_properties)
    if len(invalid_properties) > 0:
        raise FormulaError(
            # The attribute name starts with "_", so we need to remove it.
            f"`{attribute.name[1:]}` contains invalid meta properties: "
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


@frozen
class TraitFormula:
    """The trait formula.

    Attributes:
        description (str): The description of the trait.
        name (str): The name of the trait.
        type (str): The type of the trait. Either "structure" or "composition".
        coefficient (float): The coefficient of the trait.
        sia_linkage (bool): Whether the formula contains sia linkage meta properties.

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
    _numerator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_numerator],
    )
    _denominator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_denominator],
    )
    coefficient: float = field(default=1.0, validator=attrs.validators.gt(0))

    sia_linkage: bool = field(init=False, default=False)
    _initialized = field(init=False, default=False)
    _numerator = field(init=False, default=None)
    _denominator = field(init=False, default=None)

    def __attrs_post_init__(self):
        object.__setattr__(self, "sia_linkage", self._init_sia_linkage())

    def _init_sia_linkage(self) -> bool:
        """Whether the formula contains sia linkage meta properties."""
        sia_meta_properties = available_meta_properties(
            self.type, sia_linkage=True, only_sia_linkage=True  # type: ignore
        )
        for prop in itertools.chain(
            self.numerator_properties, self.denominator_properties
        ):
            if prop in sia_meta_properties:
                return True
        return False

    @property
    def numerator_properties(self) -> list[str]:
        """The meta properties in the numerator."""
        return self._numerator_properties.copy()

    @property
    def denominator_properties(self) -> list[str]:
        """The meta properties in the denominator."""
        return self._denominator_properties.copy()

    def initialize(self, meta_property_table: pd.DataFrame) -> None:
        """Initialize the trait formula.

        Args:
            meta_property_table (pd.DataFrame): The table of meta properties generated
                by `build_meta_property_table`.
        """
        numerator = self._initialize(meta_property_table, self.numerator_properties)
        object.__setattr__(self, "_numerator", numerator)
        denominator = self._initialize(meta_property_table, self.denominator_properties)
        object.__setattr__(self, "_denominator", denominator)
        object.__setattr__(self, "_initialized", True)

    @staticmethod
    def _initialize(
        meta_property_table: pd.DataFrame, properties: list[str]
    ) -> NDArray:
        return np.asarray(meta_property_table[properties].prod(axis=1))

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
        denominator[denominator == 0] = np.nan
        return numerator / denominator * self.coefficient

    def is_child_of(self, other: TraitFormula) -> bool:
        """Whether this formula is a child of the other formula.

        A formula is a child of another formula
        if both numerator and denominator of this formula have the same additional meta
        property as the other formula.

        For example,
        - A2FG is a child of A2G
        - A2FSG is a child of A2FG, and also a child of A2SG.

        Note that A2FSG is not a child of A2G, for it is a grandchild of A2G.
        Also, A2FG is not a child of CG.
        """
        if self.name in ("A1S", "A2S", "A3S", "A4S") and other.name == "CS":
            return True
        if self.name in ("A1E", "A2E", "A3E", "A4E") and other.name == "CE":
            return True
        if self.name in ("A1L", "A2L", "A3L", "A4L") and other.name == "CL":
            return True

        num1 = set(self.numerator_properties)
        den1 = set(self.denominator_properties)
        try:
            num2 = set(other.numerator_properties)
            den2 = set(other.denominator_properties)
        except AttributeError:
            raise TypeError("The other formula is not a TraitFormula instance.")

        # General case
        condition = (
            not (num2 - num1)
            and not (den2 - den1)
            and num1 - num2 == den1 - den2
            and len(num1 - num2) == 1
        )

        return condition


def load_formulas(
    type: Literal["structure", "composition"], user_file: Optional[str] = None
) -> list[TraitFormula]:
    """Load both the default formulas and the user-defined formulas.

    Args:
        type (Literal["structure", "composition"]): The type of the formulas.
        user_file (Optional[str], optional): The path of the user-defined formula file.

    Returns:
        list[TraitFormula]: The formulas.

    Raises:
        FormulaError: If a formula string cannot be parsed,
            or the user-provided formula file is in a wrong format.
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

    Raises:
        FormulaError: If a formula string cannot be parsed,
            or the user-provided formula file is in a wrong format.
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
                    numerator_properties=num_prop,  # type: ignore
                    denominator_properties=den_prop,  # type: ignore
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

    Raises:
        FormulaError: If the expression is invalid.
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


def save_trait_formula_template(dirpath: str | Path) -> None:
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


def save_builtin_formula(dirpath: str | Path) -> None:
    """Copy the builtin formula file to the given path.

    Args:
        dirpath (str): The path to save the built-in formula file.
    """
    Path(dirpath).mkdir(parents=True, exist_ok=True)
    struc_file = Path(dirpath) / "struc_builtin_formulas.txt"
    comp_file = Path(dirpath) / "comp_builtin_formulas.txt"
    struc_content = default_struc_formula_file.open("r").read()
    comp_content = default_comp_formula_file.open("r").read()
    with open(struc_file, "w", encoding="utf8") as f:
        f.write(struc_content)
    with open(comp_file, "w", encoding="utf8") as f:
        f.write(comp_content)
