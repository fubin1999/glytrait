"""Preprocess the abundance table.

Functions:
    preprocess_pipeline: Preprocess the abundance table.
    filter_glycans: Filter glycans with too many missing values.
    impute: Impute the missing values.
    normalization: Normalize the abundance table by dividing the sum of each sample.
"""
from typing import Literal

import pandas as pd


def preprocess_pipeline(
    abundance_df: pd.DataFrame,
    filter_max_na: float,
    impute_method: Literal["zero", "min", "lod", "mean", "median"],
) -> pd.DataFrame:
    """Preprocess the abundance table.

    The following steps will be performed:
        - Filter glycans with too many missing values.
        - Impute the missing values.
        - Normalize the abundance table by dividing the sum of each sample.

    Args:
        abundance_df (pd.DataFrame): The abundance table, with samples as index and Compositions
            as columns.
        filter_max_na (float): The maximum proportion of missing values allowed for a glycan.
        impute_method (str): The imputation method. Can be "zero", "min", "1/5min", "mean",
            "median".

    Returns:
        preprocessed_df (pd.DataFrame): The preprocessed abundance table.
    """
    filtered_df = filter_glycans(abundance_df, filter_max_na)
    imputed_df = impute(filtered_df, impute_method)
    preprocessed_df = normalization(imputed_df)
    return preprocessed_df


def filter_glycans(abundance_df: pd.DataFrame, max_na: float) -> pd.DataFrame:
    """Filter glycans with too many missing values.

    Args:
        abundance_df (pd.DataFrame): The abundance table, with samples as index and Compositions
            as columns.
        max_na (float): The maximum proportion of missing values allowed for a glycan.

    Returns:
        filtered_df (pd.DataFrame): The filtered abundance table.
    """
    filtered_df = abundance_df.loc[:, abundance_df.isna().mean() <= max_na]
    return filtered_df.copy()


def impute(
    abundance_df: pd.DataFrame,
    method: Literal["zero", "min", "lod", "mean", "median"],
) -> pd.DataFrame:
    """Impute the missing values.

    The following imputation methods supported:
        - "zero": Replace the missing values with 0.
        - "min": Replace the missing values with the minimum value of the corresponding glycan.
        - "lod": Replace the missing values with 1/5 of the minimum value of the corresponding
            glycan.
        - "mean": Replace the missing values with the mean value of the corresponding glycan.
        - "median": Replace the missing values with the median value of the corresponding glycan.

    Args:
        abundance_df (pd.DataFrame): The abundance table, with samples as index and Compositions
            as columns.
        method (str): The imputation method. Can be "zero", "min", "1/5min", "mean", "median".

    Returns:
        imputed_df (pd.DataFrame): The imputed abundance table.
    """
    if method == "zero":
        imputed_df = abundance_df.fillna(0)
    elif method == "min":
        imputed_df = abundance_df.fillna(abundance_df.min())
    elif method == "lod":
        imputed_df = abundance_df.fillna(abundance_df.min() / 5)
    elif method == "mean":
        imputed_df = abundance_df.fillna(abundance_df.mean())
    elif method == "median":
        imputed_df = abundance_df.fillna(abundance_df.median())
    else:
        raise ValueError("The imputation method is not supported.")
    return imputed_df


def normalization(abundance_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize the abundance table by dividing the sum of each sample.

    Args:
        abundance_df (pd.DataFrame): The abundance table, with samples as index and Compositions
            as columns.

    Returns:
        normalized_df (pd.DataFrame): The normalized abundance table.
    """
    row_sums = abundance_df.sum(axis=1)
    normalized_df = abundance_df.div(row_sums, axis=0).round(6)
    return normalized_df
