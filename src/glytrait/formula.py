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
from collections.abc import Callable
from importlib.resources import files, as_file
from pathlib import Path
from typing import Literal, Optional, Generator, TypeVar, Type

import numpy as np
import pandas as pd
import attrs
from attrs import field, define

from glytrait.data_type import AbundanceTable, MetaPropertyTable

__all__ = [
    "TraitFormula",
    "load_formulas",
    "save_builtin_formula",
]

from glytrait.exception import GlyTraitError

default_struc_formula_file = files("glytrait.resources").joinpath("struc_formula.txt")
default_comp_formula_file = files("glytrait.resources").joinpath("comp_formula.txt")


class FormulaError(GlyTraitError):
    """Raised if a formula is invalid."""


class FormulaFileError(GlyTraitError):
    """Raised if a formula file is in wrong format."""


class FormulaParseError(FormulaError):
    """Raised when failing parsing a formula expression."""

    def __init__(self, expr: str, reason: str = ""):
        super().__init__()
        self.expr = expr
        self.reason = reason

    def __str__(self):
        return f"Invalid expression: {self.expr}. {self.reason}"


class FormulaTermParseError(FormulaParseError):
    """Raised when failing parsing a formula term expression."""


class FormulaCalculationError(FormulaError):
    """Raised when failing to calculate a derived trait."""


class FormulaTermCalculationError(FormulaCalculationError):
    """Raised when calling a formula term fails."""

    def __init__(self, term: str):
        super().__init__()
        self.term = term


class MissingMetaPropertyError(FormulaTermCalculationError):
    """Raised when a meta-property is missing in the meta-property table."""

    def __init__(self, term: str, mp: str):
        super().__init__(term)
        self.mp = mp

    def __str__(self):
        return f"Meta-property '{self.mp}' is missing in the meta-property table."


class MetaPropertyTypeError(FormulaTermCalculationError):
    """Raised when the dtype of a meta-property is not supported."""

    def __init__(self, term: str, mp: str, dtype: str, reason: str = ""):
        super().__init__(term)
        self.mp = mp
        self.dtype = dtype
        self.reason = reason

    def __str__(self):
        return (
            f"Meta-property '{self.mp}' has an unsupported dtype: "
            f"'{self.dtype}'. {self.reason}"
        )


class FormulaNotInitializedError(FormulaError):
    """Raised when a formula is not initialized."""


TTerm = TypeVar("TTerm", bound="FormulaTerm")


class FormulaTerm:
    """The protocol for a formula term.

    A formula term is a callable object that takes a meta-property table as input
    and returns a Series with the same index as the meta-property table.
    It makes some calculations based on the meta-properties.
    Normally, it will only work with one meta property.

    A common using paradigm is to:
    1. Use `meets_condition` to check if the expression could be parsed into the term.
    2. Use `from_expr` to create a formula term from an expression.
    """

    def __call__(self, meta_property_table: MetaPropertyTable) -> pd.Series:
        """Calculate the term.

        The return value is a Series with the same index as the meta-property table.
        (The index is the glycans.)

        Args:
            meta_property_table: The table of meta properties.

        Returns:
            pd.Series: The values of the term,
                with the same index as the meta-property table.

        Raises:
            FormulaCalculationError: If calling the term fails.
        """
        raise NotImplementedError

    @classmethod
    def from_expr(cls: Type[TTerm], expr: str) -> TTerm:
        """Create a formula term from an expression.

        Args:
            expr: The expression of the term.

        Returns:
            The formula term.

        Raises:
            FormulaTermParseError: If the expression is invalid.
        """
        raise NotImplementedError

    @staticmethod
    def meets_condition(expr: str) -> bool:
        """Check if the expression could be parsed into the term.

        Args:
            expr: The condition expression.

        Returns:
            Whether the term meets the condition.
        """
        raise NotImplementedError

    @property
    def expr(self) -> str:
        """The expression of the term."""
        raise NotImplementedError


_terms: list[type[FormulaTerm]] = []


def _register_term(cls: type[FormulaTerm]) -> type[FormulaTerm]:
    """A decorator to register a formula term class."""
    _terms.append(cls)
    return cls


@_register_term
@define
class ConstantTerm(FormulaTerm):
    """Return a series with all values being constant."""

    value: int = field()

    def __call__(self, meta_property_table: MetaPropertyTable) -> pd.Series:
        """Calculate the term."""
        # Add validation for the dtype of column
        return pd.Series(
            self.value, index=meta_property_table.index, name=self.expr, dtype="UInt8"
        )

    @property
    def expr(self) -> str:
        """The expression of the term."""
        return str(self.value)

    @staticmethod
    def meets_condition(expr: str) -> bool:
        """Check if the expression could be parsed into the term."""
        return expr.strip("()").isdigit()

    @classmethod
    def from_expr(cls, expr: str) -> ConstantTerm:
        """Create a formula term from an expression."""
        if not cls.meets_condition(expr):
            raise FormulaTermParseError(expr, "Conditions not met.")
        return cls(value=int(expr.strip("()")))


@_register_term
@define
class NumericalTerm(FormulaTerm):
    """Return the values of a numerical meta-property.

    This term simply returns a column of the meta-property table,
    as a Series with the same index as the meta-property table.

    Note that the dtype of the meta-property should be UInt8.

    Args:
        meta_property: The numerical meta property.
    """

    meta_property: str = field()

    def __call__(self, meta_property_table: MetaPropertyTable) -> pd.Series:
        """Calculate the term.

        The return value is a Series with the same index as the meta-property table.

        Args:
            meta_property_table: The table of meta properties.

        Returns:
            pd.Series: The values of the meta-property,
            with the same index as the meta-property table.

        Raises:
            MissingMetaPropertyError: If the meta-property is not in the meta-property table.
            MetaPropertyTypeError: If the dtype of the meta-property is not UInt8.
        """
        try:
            mp_s = meta_property_table[self.meta_property]
        except KeyError as e:
            raise MissingMetaPropertyError(self.expr, self.meta_property) from e

        if mp_s.dtype == "boolean" or mp_s.dtype == "category":
            reason = "NumericalTerm only works with numerical meta properties."
            raise MetaPropertyTypeError(self.meta_property, str(mp_s.dtype), reason)

        return mp_s

    @property
    def expr(self) -> str:
        """The expression of the term."""
        return self.meta_property

    @staticmethod
    def meets_condition(expr: str) -> bool:
        """Check if the expression could be parsed into the term."""
        # Consist of digits, letters and "_", but not all digits
        pattern = r"\w+"
        return bool(re.fullmatch(pattern, expr.strip("()"))) and not expr.isdigit()

    @classmethod
    def from_expr(cls, expr: str) -> NumericalTerm:
        """Create a formula term from an expression."""
        if not cls.meets_condition(expr):
            raise FormulaTermParseError(expr, "Conditions not met.")
        return cls(meta_property=expr.strip("()"))


OPERATORS = {"==", "!=", ">", ">=", "<", "<="}


@_register_term
@define
class CompareTerm(FormulaTerm):
    """Compare the value of a meta-property with a given value.

    The validity of `meta_property` will not be checked here.

    Args:
        meta_property: The meta property to compare.
        operator: The comparison operator.
        value: The value to compare with.
    """

    meta_property: str = field()
    operator: Literal["==", "!=", ">", ">=", "<", "<="] = field()
    value: float | bool | str = field()

    def __call__(self, meta_property_table: MetaPropertyTable) -> pd.Series:
        """Calculate the term.

        Args:
            meta_property_table: The table of meta properties.

        Returns:
            pd.Series: A boolean Series with the same index as the meta-property table.

        Raises:
            MissingMetaPropertyError: If the meta-property is not in the meta-property table.
            MetaPropertyTypeError: If the dtype of the meta-property is not supported.
        """
        try:
            mp_s = meta_property_table[self.meta_property]
        except KeyError as e:
            raise MissingMetaPropertyError(self.expr, self.meta_property) from e

        condition_1 = mp_s.dtype == "boolean" or mp_s.dtype == "category"
        condition_2 = self.operator in {">", ">=", "<", "<="}
        if condition_1 and condition_2:
            reason = (
                f"Cannot use '{self.operator}' with '{mp_s.dtype}' meta properties."
            )
            raise MetaPropertyTypeError(
                self.expr, self.meta_property, str(mp_s.dtype), reason
            )

        if isinstance(self.value, str):
            expr = f"mp_s {self.operator} '{self.value}'"
        else:
            expr = f"mp_s {self.operator} {self.value}"
        result_s = eval(expr)

        result_s.name = self.expr
        result_s = result_s.astype("UInt8")
        return result_s

    @property
    def expr(self) -> str:
        """The expression of the term."""
        if isinstance(self.value, str):
            return f"{self.meta_property} {self.operator} '{self.value}'"
        else:
            return f"{self.meta_property} {self.operator} {self.value}"

    @staticmethod
    def meets_condition(expr: str) -> bool:
        """Check if the expression could be parsed into the term."""
        for op in OPERATORS:
            if op in expr:
                return True
        return False

    @classmethod
    def from_expr(cls, expr: str) -> CompareTerm:
        """Create a formula term from an expression."""
        if not cls.meets_condition(expr):
            raise FormulaTermParseError(expr, "Conditions not met.")
        if not expr.startswith("(") or not expr.endswith(")"):
            raise FormulaTermParseError(expr, "Missing parentheses.")
        expr = expr.strip("()")

        meta_property_p = r"(\w+)"
        operator_p = r"(==|!=|>|>=|<|<=)"
        value_p = r"""(\d+|True|False|'[a-zA-Z-_]*'|"[a-zA-Z-_]*")"""
        total_p = rf"{meta_property_p}\s*{operator_p}\s*{value_p}"
        match = re.fullmatch(total_p, expr)
        if match is None:
            raise FormulaTermParseError(expr, "Unknown error.")
        meta_property, operator, value = match.groups()

        if value.isdigit():
            value = int(value)
        elif value == "True":
            value = True
        elif value == "False":
            value = False
        else:
            value = value.strip("'")
            value = value.strip('"')

        return cls(  # type: ignore
            meta_property=meta_property, operator=operator, value=value  # type: ignore
        )


def _is_sia_linkage_term(term: FormulaTerm) -> bool:
    """Check if a formula term is related to sia-linkage."""
    return "nE" in term.expr or "nL" in term.expr


FormulaParserType = Callable[[str], tuple[str, list[FormulaTerm], list[FormulaTerm]]]


@define
class TraitFormula:
    """The trait formula."""

    name: str = field()
    description: str = field()
    numerators: list[FormulaTerm] = field(converter=list, validator=attrs.validators.min_len(1))
    denominators: list[FormulaTerm] = field(converter=list, validator=attrs.validators.min_len(1))

    _sia_linkage: bool = field(init=False)
    _initialized: bool = field(init=False, default=False)
    _numerator_array: pd.Series = field(init=False, default=None)
    _denominator_array: pd.Series = field(init=False, default=None)

    def __attrs_post_init__(self):
        self._sia_linkage = self._init_sia_linkage()

    def _init_sia_linkage(self) -> bool:
        """Check if the formula is related to sia-linkage."""
        term_iter = itertools.chain(self.numerators, self.denominators)
        return any(_is_sia_linkage_term(term) for term in term_iter)

    @classmethod
    def from_expr(
        cls,
        expr: str,
        description: str,
        *,
        parser: Optional[FormulaParserType] = None  # Dependency injection
    ) -> TraitFormula:
        """Create a formula from an expression.

        Args:
            expr: The expression of the formula.
            description: The description of the formula.
            parser: The parser to parse the formula expression.
                If None, the default parser will be used.
                Defaults to None.

        Returns:
            The formula.

        Raises:
            FormulaParseError: If the expression is invalid.
        """
        if parser is None:
            parser = _parse_formula_expression
        name, numerators, denominators = parser(expr)
        return cls(
            name=name,
            description=description,
            numerators=numerators,
            denominators=denominators
        )

    @property
    def sia_linkage(self) -> bool:
        """Whether the formula is related to sia-linkage."""
        return self._sia_linkage

    def initialize(self, meta_property_table: MetaPropertyTable) -> None:
        """Initialize the trait formula.

        Args:
            meta_property_table: The table of meta properties.

        Raises:
            FormulaCalculationError: If the initialization fails.
        """
        self._numerator_array = self._initialize(meta_property_table, self.numerators)
        self._denominator_array = self._initialize(
            meta_property_table, self.denominators
        )
        self._initialized = True

    @staticmethod
    def _initialize(
        meta_property_table: MetaPropertyTable, terms: list[FormulaTerm]
    ) -> pd.Series:
        """Initialize the numerator or denominator of the formula."""
        try:
            series_list = [term(meta_property_table) for term in terms]
        except FormulaTermCalculationError as e:
            msg = f"Failed to calculate term: {e.term}. {str(e)}"
            raise FormulaCalculationError(msg) from e
        return pd.concat(series_list, axis=1).prod(axis=1)

    def calcu_trait(self, abundance_table: AbundanceTable) -> pd.Series:
        """Calculate the trait.

        Args:
            abundance_table (AbundanceTable): The glycan abundance table,
                with samples as index, and glycans as columns.

        Returns:
            pd.Series: An array of trait values for each sample.

        Raises:
            FormulaNotInitializedError: If the formula is not initialized.
        """
        if not self._initialized:
            raise FormulaNotInitializedError()
        pd.testing.assert_index_equal(
            abundance_table.columns, self._numerator_array.index
        )
        pd.testing.assert_index_equal(
            abundance_table.columns, self._denominator_array.index
        )

        numerator = abundance_table.values @ self._numerator_array.values
        denominator = abundance_table.values @ self._denominator_array.values
        denominator[denominator == 0] = np.nan
        values = numerator / denominator
        return pd.Series(
            values, index=abundance_table.index, name=self.name, dtype=float
        )


def _parse_formula_expression(
    expr: str,
) -> tuple[str, list[FormulaTerm], list[FormulaTerm]]:
    """Parse the formula expression.

    Args:
        expr: The formula expression.

    Returns:
        The name, the numerators (a list of FormulaTerm),
        and the denominators (a list of FormulaTerm) of the formula.

    Raises:
        FormulaParseError: If the expression is invalid.
    """
    name, numerator_expr, spliter, denominator_expr = _split_formula_expression(expr)
    numerators = _parse_terms(numerator_expr)
    denominators = _parse_terms(denominator_expr)
    if spliter == "//":
        numerators.extend(denominators)
    return name, numerators, denominators


def _split_formula_expression(expr: str) -> tuple[str, str, str, str]:
    """Split the formula expression into four parts:

    - The name of the formula.
    - The numerator expression of the formula.
    - The slash or double-slash.
    - The denominator expression of the formula.

    Args:
        expr: The formula expression.

    Returns:
        The name, the numerator expression, the spliter, and the denominator expression.

    Raises:
        FormulaParseError: If the expression is invalid.
    """
    expr = expr.strip()
    if " = " not in expr:
        raise FormulaParseError(expr, "Missing ' = ' (spaces are important).")

    try:
        name, expr_after_name = expr.split(" = ")
    except ValueError:
        raise FormulaParseError(expr, "Misuse of '=' for '=='.")
    name = name.strip()

    if "//" not in expr_after_name and "/" not in expr_after_name:
        raise FormulaParseError(expr, "no '/' or '//'.")

    if expr_after_name.count("//") == 1:
        spliter = "//"
    elif expr_after_name.count("/") == 1:
        spliter = "/"
    else:
        raise FormulaParseError(expr, "too many '/' or '//'.")
    numerator, denominator = expr_after_name.split(spliter)
    numerator = numerator.strip()
    denominator = denominator.strip()

    return name, numerator, spliter, denominator


def _parse_terms(expr: str) -> list[FormulaTerm]:
    """Parse the terms in the formula expression.

    Args:
        expr: The expression of the terms.

    Returns:
        The terms.

    Raises:
        FormulaParseError: If the expression is invalid.
    """
    term_exprs = [t.strip() for t in expr.split("*")]
    try:
        return [_parse_term(t) for t in term_exprs]
    except FormulaTermParseError as e:
        reason = f"Could not parse term: {e.expr}. {e.reason}"
        raise FormulaParseError(expr, reason) from e


def _parse_term(expr: str) -> FormulaTerm:
    """Parse a term in the formula expression.

    Args:
        expr: The expression of the term.

    Returns:
        The formula term.

    Raises:
        FormulaTermParseError: If the expression is invalid.
    """
    for term_cls in _terms:
        if term_cls.meets_condition(expr):
            return term_cls.from_expr(expr)  # Let the exception pass through
    raise FormulaTermParseError(expr, "Does not belong to any term class.")


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
        FormulaFileError: If the formula file is in a wrong format.
        FormulaParseError: If a formula string cannot be parsed.
    """
    formulas = list(load_default_formulas(type_=type_))

    if user_file is not None:
        default_formula_names = {f.name for f in formulas}
        user_formulas = load_formulas_from_file(user_file)
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
        yield from load_formulas_from_file(str(file))


def load_formulas_from_file(
    filepath: str,
    *,
    parser: Optional[FormulaParserType] = None,  # Dependency injection
) -> Generator[TraitFormula, None, None]:
    """Load formulas from a file.

    Args:
        filepath: The path of the formula file.
        parser: The parser to parse the formula expression.
            If None, the default parser will be used.
            Defaults to None.

    Yields:
        The formulas parsed.

    Raises:
        FormulaFileError: If the formula file is in a wrong format.
        FormulaParseError: If a formula string cannot be parsed.
    """
    formulas_parsed: set[str] = set()
    for description, expression in deconvolute_formula_file(filepath):
        formula = TraitFormula.from_expr(expression, description, parser=parser)
        if formula.name in formulas_parsed:
            raise FormulaFileError(f"Duplicate formula name: {formula.name}.")
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
        FormulaFileError: If the formula file is in a wrong format.
    """
    description = None
    expression = None
    with open(formula_file, "r", encoding="utf8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                if description is not None:
                    raise FormulaFileError(
                        f"No expression follows description '{description}'."
                    )
                description = line[1:].strip()
            elif line.startswith("$"):
                expression = line[1:].strip()
                if description is None:
                    raise FormulaFileError(
                        f"No description before expression '{expression}'."
                    )
                yield description, expression
                description = None
                expression = None
    if description is not None:
        raise FormulaFileError(f"No expression follows description '{description}'.")


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
