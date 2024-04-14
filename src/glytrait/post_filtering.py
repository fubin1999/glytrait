"""Post-filtering the derived traits table.

The only function in this module is `post_filter`.
"""

from typing import Iterable, Literal

import numpy as np

from glytrait.formula import TraitFormula
from glytrait.data_type import DerivedTraitTable

__all__ = ["post_filter"]


def post_filter(
    formulas: Iterable[TraitFormula],
    trait_df: DerivedTraitTable,
    threshold: float,
    method: Literal["pearson", "spearman"],
) -> DerivedTraitTable:
    """Post-filter the derived traits table.

    The post-filtering consists of two steps:
    1. Rule out the invalid traits.
    2. Reduce the colinearity of the traits.

    Args:
        formulas (Iterable[TraitFormula]): The formulas to be filtered.
        trait_df (DerivedTraitTable): The derived traits table.
        threshold (float): The threshold of the correlation coefficient.
            If set to -1, the colinearity filtering will be skipped.
        method (Literal["pearson", "spearman"]): The method to calculate the correlation.

    Returns:
        DerivedTraitTable: The filtered derived traits table.
    """
    trait_df = filter_invalid(trait_df)
    if threshold != -1:
        trait_df = filter_colinearity(formulas, trait_df, threshold, method)
    return trait_df


def filter_invalid(trait_df: DerivedTraitTable) -> DerivedTraitTable:
    """Rule out the invalid traits.

    A trait is invalid if it:
    1. Has the same value for all samples.
    2. Is NaN for all samples.

    Args:
        trait_df (DerivedTraitTable): The derived traits table.

    Returns:
        DerivedTraitTable: The filtered trait values.
    """
    trait_df = _filter_all_same(trait_df)
    trait_df = _filter_all_nan(trait_df)
    return trait_df


def _filter_all_same(trait_df: DerivedTraitTable) -> DerivedTraitTable:
    """Rule out the traits that have the same value for all samples."""
    return DerivedTraitTable(trait_df.loc[:, trait_df.nunique() != 1])


def _filter_all_nan(trait_df: DerivedTraitTable) -> DerivedTraitTable:
    """Rule out the traits that are NaN for all samples."""
    return DerivedTraitTable(trait_df.loc[:, trait_df.notna().any()])


def filter_colinearity(
    formulas: Iterable[TraitFormula],
    trait_df: DerivedTraitTable,
    threshold: float,
    method: Literal["pearson", "spearman"],
) -> DerivedTraitTable:
    """Filter the colinearity of the formulas.

    Args:
        formulas (Iterable[TraitFormula]): The formulas to be filtered.
        trait_df (DerivedTraitTable): The derived traits table, after post-filtering.
        threshold (float): The threshold of the correlation coefficient.
        method (Literal["pearson", "spearman"]): The method to calculate the correlation.

    Returns:
        DerivedTraitTable: The filtered derived traits table.
    """
    trait_names = list(trait_df.columns)

    # First, build the parent-child relationship matrix.
    # This matrix consists of 0 and 1.
    # If the i-th row and j-th column is 1, then the i-th trait is the child of the j-th trait.
    rela_matrix = _relationship_matrix(trait_names, formulas)

    # Second, build the correlation matrix.
    # This matrix also consists of 0 and 1.
    # 1 means the two traits are highly correlated (r > the threshold).
    corr_matrix = _correlation_matrix(trait_df, threshold, method)

    # Third, multiply the two matrices.
    # If the i-th row and j-th column is 1,
    # then the i-th trait is the child of the j-th trait,
    # and the two traits are highly correlated.
    # Then the i-th trait should be removed.
    # That is to say, if the row sum of the i-th row is larger than 0,
    # then the i-th trait should be removed.
    remove_matrix = rela_matrix * corr_matrix
    to_keep = np.sum(remove_matrix, axis=1) == 0

    # Finally, filter the formulas and the derived traits table.
    filtered_trait_table = trait_df.loc[:, to_keep]
    return DerivedTraitTable(filtered_trait_table)


def _relationship_matrix(trait_names: Iterable[str], formulas: Iterable[TraitFormula]):
    """Build the parent-child relationship matrix.

    This matrix consists of 0 and 1.
    If the i-th row and j-th column is 1, then the i-th trait is the child of the j-th trait.

    Args:
        trait_names (Iterable[str]): The names of the traits.
        formulas (Iterable[TraitFormula]): The formulas to be filtered.

    Returns:
        np.ndarray: The relationship matrix.
    """
    trait_names = list(trait_names)
    formula_dict = {
        formula.name: formula for formula in formulas if formula.name in trait_names
    }
    matrix = np.zeros((len(trait_names), len(trait_names)), dtype=int)
    for i, trait_1 in enumerate(trait_names):
        for j, trait_2 in enumerate(trait_names):
            formula_1 = formula_dict[trait_1]
            formula_2 = formula_dict[trait_2]
            if _is_child_of(formula_1, formula_2):
                matrix[i, j] = 1
    return matrix


def _correlation_matrix(
    trait_table: DerivedTraitTable,
    threshold: float,
    method: Literal["pearson", "spearman"],
):
    """Build the correlation matrix.

    This matrix also consists of 0 and 1.
    1 means the two traits are highly correlated (r > the threshold).

    Args:
        trait_table (DerivedTraitTable): The derived traits table, after post-filtering.
        threshold (float): The threshold of the correlation coefficient.
        method (Literal["pearson", "spearman"]): The method to calculate the correlation.

    Returns:
        np.ndarray: The correlation matrix.
    """
    corr_matrix = trait_table.corr(method=method).values >= threshold
    corr_matrix = corr_matrix.astype(int)
    return corr_matrix


def _is_child_of(trait1: TraitFormula, trait2: TraitFormula) -> bool:
    """Whether trait1 is a child of the trait2.

    A formula is a child of another formula
    if both numerator and denominator of this formula have the same additional
    term as the other formula.

    For example,
    Formula 1:
      - Numerators: A, B, C
      - Denominators: A, B
    Formula 2:
      - Numerators: A, B, C, D
      - Denominators: A, B, D
    Then Formula 2 is a child of Formula 1.

    """
    num1 = set(t.expr for t in trait1.numerators)
    den1 = set(t.expr for t in trait1.denominators)
    try:
        num2 = set(t.expr for t in trait2.numerators)
        den2 = set(t.expr for t in trait2.denominators)
    except AttributeError:
        raise TypeError("The other formula is not a TraitFormula instance.")

    condition = (
        not (num2 - num1)
        and not (den2 - den1)
        and num1 - num2 == den1 - den2
        and len(num1 - num2) == 1
    )

    return condition
