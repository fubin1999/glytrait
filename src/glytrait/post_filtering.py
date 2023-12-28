from typing import Iterable, Literal

import numpy as np

from glytrait.formula import TraitFormula
from glytrait.trait import DerivedTraitTable


def filter_invalid(
    formulas: Iterable[TraitFormula], trait_df: DerivedTraitTable
) -> tuple[list[TraitFormula], DerivedTraitTable]:
    """Rule out the invalid traits.

    A trait is invalid if it:
    1. Has the same value for all samples.
    2. Is NaN for all samples.

    Args:
        formulas (Iterable[TraitFormula]): The formulas to be filtered.
        trait_df (DerivedTraitTable): The trait values, with samples as index and trait names
            as columns.

    Returns:
        DerivedTraitTable: The filtered trait values.
    """
    trait_df = _filter_all_same(trait_df)
    trait_df = _filter_all_nan(trait_df)
    formulas = [f for f in formulas if f.name in trait_df.columns]
    return formulas, trait_df


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
) -> tuple[list[TraitFormula], DerivedTraitTable]:
    """Filter the colinearity of the formulas.

    Args:
        formulas (Iterable[TraitFormula]): The formulas to be filtered.
        trait_df (DerivedTraitTable): The derived traits table, after post-filtering.
        threshold (float): The threshold of the correlation coefficient.
        method (Literal["pearson", "spearman"]): The method to calculate the correlation.

    Returns:
        tuple[list[TraitFormula],DerivedTraitTable]: The filtered formulas and the filtered
            derived traits table.
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
    formulas = [f for f in formulas if f.name in filtered_trait_table.columns]
    return formulas, DerivedTraitTable(filtered_trait_table)


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
