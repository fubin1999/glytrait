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
from importlib.resources import files, as_file
from pathlib import Path
from typing import Literal, Optional, Generator

import numpy as np
import pandas as pd
from attrs import field, frozen, validators

from glytrait.exception import FormulaError
from glytrait.data_type import AbundanceTable, MetaPropertyTable
from glytrait.meta_property import available_meta_properties

__all__ = [
    "TraitFormula",
    "load_formulas",
    "save_trait_formula_template",
    "save_builtin_formula",
]

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

    The `TraitFormula` is used to represent a trait formula,
    and to calculate the trait values for each sample.
    To use it, you need to
    1. Make sure the index of the meta-property table (glycans) and the columns of the
    abundance table (also glycans) are in the same order.
    2. Initialize the formula by calling `initialize` method with a meta-property table.
    3. Calculate the trait values by calling `calcu_trait` method with an abundance table.

    Notes:
        Please make sure that the index of the meta-property table
        and the columns of the abundance table are in the same order.
        Otherwise, an assertion error will be raised.
        This can be easily achieved by using the `reindex` method of pandas DataFrame
        (see the example below).
        The reason why `TraitFormula` does not do this automatically is for efficiency.
        Reindexing the meta-property table on each call of `calcu_trait` is time-consuming
        and unnecessary.

    Attributes:
        description (str): The description of the trait.
        name (str): The name of the trait.
        type (str): The type of the trait. Either "structure" or "composition".
        numerator_properties (list[str]): The meta properties in the numerator.
        denominator_properties (list[str]): The meta properties in the denominator.
        coefficient (float): The coefficient of the trait.
        sia_linkage (bool): Whether the formula contains sia linkage meta properties.

    Examples:
        >>> from glytrait.formula import TraitFormula
        >>> meta_property_table = ...
        >>> abundance_table = ...
        >>> formula = TraitFormula(
        ...     description="The ratio of high-mannose to complex glycans",
        ...     name="MHy",
        ...     type="structure",
        ...     numerator_properties=["isHighMannose"],
        ...     denominator_properties=["isComplex"],
        ... )
        >>> meta_property_table = meta_property_table.reindex(abundance_table.columns)
        >>> pd.testing.assert_index_equal(meta_property_table.index, abundance_table.columns)
        >>> formula.initialize(meta_property_table)
        >>> trait_values = formula.calcu_trait(abundance_table)  # a Series
    """

    description: str = field()
    name: str = field()
    type: str = field(validator=validators.in_(["structure", "composition"]))
    _numerator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_numerator],
    )
    _denominator_properties: list[str] = field(
        converter=list,
        validator=[_check_length, _check_meta_properties, _check_denominator],
    )
    coefficient: float = field(default=1.0, validator=validators.gt(0))

    sia_linkage: bool = field(init=False, default=False)
    _initialized: bool = field(init=False, default=False)
    _numerator: pd.Series = field(init=False, default=None)
    _denominator: pd.Series = field(init=False, default=None)

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

    def initialize(self, meta_property_table: MetaPropertyTable) -> None:
        """Initialize the trait formula.

        Args:
            meta_property_table (MetaPropertyTable):
                The table of meta properties generated by `build_meta_property_table`.
        """
        numerator = self._initialize(meta_property_table, self.numerator_properties)
        object.__setattr__(self, "_numerator", numerator)
        denominator = self._initialize(meta_property_table, self.denominator_properties)
        object.__setattr__(self, "_denominator", denominator)
        object.__setattr__(self, "_initialized", True)

    @staticmethod
    def _initialize(
        meta_property_table: MetaPropertyTable, properties: list[str]
    ) -> pd.Series:
        return meta_property_table[properties].prod(axis=1)

    def calcu_trait(self, abundance_table: AbundanceTable) -> pd.Series:
        """Calculate the trait.

        Args:
            abundance_table (AbundanceTable): The glycan abundance table,
                with samples as index, and glycans as columns.

        Returns:
            pd.Series: An array of trait values for each sample.
        """
        if not self._initialized:
            raise RuntimeError("TraitFormula is not initialized.")
        pd.testing.assert_index_equal(abundance_table.columns, self._numerator.index)
        pd.testing.assert_index_equal(abundance_table.columns, self._denominator.index)

        numerator = abundance_table.values @ self._numerator.values
        denominator = abundance_table.values @ self._denominator.values
        denominator[denominator == 0] = np.nan
        values = numerator / denominator * self.coefficient
        return pd.Series(
            values, index=abundance_table.index, name=self.name, dtype=float
        )

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
    type_: Literal["structure", "composition"],
    user_file: Optional[str] = None,
    sia_linkage: bool = False,
) -> list[TraitFormula]:
    """Load both the default formulas and the user-defined formulas.

    Args:
        type_ (Literal["structure", "composition"]): The type of the formulas.
        user_file (Optional[str], optional): The path of the user-defined formula file.
            Defaults to None.
        sia_linkage (bool, optional): Whether to include formulas with sia-linkage
            meta-properties. Defaults to False.

    Returns:
        list[TraitFormula]: The formulas.

    Raises:
        FormulaError: If a formula string cannot be parsed,
            or the user-provided formula file is in a wrong format.
    """
    formulas = list(load_default_formulas(type_=type_))

    if user_file is not None:
        default_formula_names = {f.name for f in formulas}
        user_formulas = load_formulas_from_file(user_file, type_=type_)
        for f in user_formulas:
            if f.name in default_formula_names:
                continue
            formulas.append(f)

    if not sia_linkage:
        formulas = [f for f in formulas if not f.sia_linkage]
    return formulas


def load_default_formulas(
    type_: Literal["structure", "composition"]
) -> Generator[TraitFormula, None, None]:
    """Load the default formulas.

    Args:
        type_ (Literal["structure", "composition"]): The type of the formulas.

    Yields:
        The formulas parsed.
    """
    if type_ == "composition":
        file_traversable = default_comp_formula_file
    elif type_ == "structure":
        file_traversable = default_struc_formula_file
    else:
        raise ValueError("Invalid formula type.")
    with as_file(file_traversable) as file:
        yield from load_formulas_from_file(str(file), type_=type_)


def load_formulas_from_file(
    filepath: str, type_: Literal["structure", "composition"]
) -> Generator[TraitFormula, None, None]:
    """Load formulas from a file.

    Args:
        filepath (str): The path of the formula file.
        type_ (Literal["structure", "composition"]): The type of the formulas.

    Yields:
        The formulas parsed.

    Raises:
        FormulaError: If a formula string cannot be parsed,
            or the user-provided formula file is in a wrong format,
            or there are duplicate formula names.
    """
    formulas_parsed: set[str] = set()
    for description, expression in deconvolute_formula_file(filepath):
        formula = create_formula(description, expression, type_=type_)
        if formula.name in formulas_parsed:
            raise FormulaError(f"Duplicate formula name: {formula.name}.")
        yield formula
        formulas_parsed.add(formula.name)


def deconvolute_formula_file(
    formula_file: str,
) -> Generator[tuple[str, str], None, None]:
    """A generator that yields the formula description and the formula expression.

    Args:
        formula_file (str): The path of the formula file.

    Yields:
        tuple[str, str]: The formula description and the formula expression.

    Raises:
        FormulaError: If the user-provided formula file is in a wrong format.
    """
    description = None
    expression = None
    with open(formula_file, "r", encoding="utf8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                if description is not None:
                    raise FormulaError(
                        f"No expression follows description '{description}'."
                    )
                description = line[1:].strip()
            elif line.startswith("$"):
                expression = line[1:].strip()
                if description is None:
                    raise FormulaError(
                        f"No description before expression '{expression}'."
                    )
                yield description, expression
                description = None
                expression = None
    if description is not None:
        raise FormulaError(f"No expression follows description '{description}'.")


def create_formula(
    description: str, expression: str, type_: Literal["structure", "composition"]
) -> TraitFormula:
    """Create a formula from the description and the expression.

    Args:
        description (str): The description of the formula.
        expression (str): The expression of the formula.
        type_ (Literal["structure", "composition"]): The type of the formula.

    Returns:
        TraitFormula: The formula.

    Raises:
        FormulaError: If the expression is invalid, or the meta properties are invalid.
    """
    name, num_prop, den_prop, coef = parse_formula_expression(expression)
    return TraitFormula(
        description=description,
        name=name,
        type=type_,
        numerator_properties=num_prop,  # type: ignore
        denominator_properties=den_prop,  # type: ignore
        coefficient=coef,
    )


def parse_formula_expression(expr: str) -> tuple[str, list[str], list[str], float]:
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
