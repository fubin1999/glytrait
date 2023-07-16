from collections.abc import Iterable

import numpy as np
import pandas as pd

from glytrait.formula import TraitFormula


def filter_colinearity(
    formulas: Iterable[TraitFormula],
    derived_traits_table: pd.DataFrame,
    threshold: float,
) -> set[TraitFormula]:
    """Filter the colinearity of the formulas.

    Args:
        formulas (Iterable[TraitFormula]): The formulas to be filtered.
        derived_traits_table (pd.DataFrame): The derived traits table, after post-filtering.
        threshold (float): The threshold of the correlation coefficient.

    Returns:
        An array of bool, indicating whether the corresponding trait should be kept.
        The array is in the order of the column of `derived_traits_table`.
    """
    trait_names = list(derived_traits_table.columns)

    # First, build the parent-child relationship matrix.
    # This matrix consists of 0 and 1.
    # If the i-th row and j-th column is 1, then the i-th trait is the child of the j-th trait.
    rela_matrix = _relationship_matrix(trait_names, formulas)

    # Second, build the correlation matrix.
    # This matrix also consists of 0 and 1.
    # 1 means the two traits are highly correlated (r > the threshold).
    corr_matrix = _correlation_matrix(derived_traits_table, threshold)

    # Third, multiply the two matrices.
    # If the i-th row and j-th column is 1,
    # then the i-th trait is the child of the j-th trait,
    # and the two traits are highly correlated.
    # Then the i-th trait should be removed.
    # That is to say, if the row sum of the i-th row is larger than 0,
    # then the i-th trait should be removed.
    remove_matrix = rela_matrix * corr_matrix
    to_keep = np.sum(remove_matrix, axis=1) == 0

    return to_keep


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
            if formula_1.is_child_of(formula_2):
                matrix[i, j] = 1
    return matrix


def _correlation_matrix(trait_table: pd.DataFrame, threshold: float):
    """Build the correlation matrix.

    This matrix also consists of 0 and 1.
    1 means the two traits are highly correlated (r > the threshold).

    Args:
        trait_table (pd.DataFrame): The derived traits table, after post-filtering.
        threshold (float): The threshold of the correlation coefficient.

    Returns:
        np.ndarray: The correlation matrix.
    """
    corr_matrix = trait_table.corr().values > threshold
    corr_matrix = corr_matrix.astype(int)
    return corr_matrix
