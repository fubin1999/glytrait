"""Preprocess the abundance table.

The only high-level function that client code will use is `preprocess`.
"""

from typing import Literal, Protocol

from attrs import define, field
from attrs.validators import in_

from glytrait.data_type import AbundanceTable

__all__ = ["preprocess"]


# This is the high-level function that client code will use.
def preprocess(
    data: AbundanceTable,
    filter_max_na: float,
    impute_method: Literal["zero", "min", "lod", "mean", "median"],
) -> AbundanceTable:
    """Preprocess the abundance table.

    Notes:
        This function will modify the input data in place.

    Args:
        data (AbundanceTable): The abundance table.
        filter_max_na (float): The maximum proportion of missing values allowed for a glycan.
        impute_method (Literal["zero", "min", "lod", "mean", "median"]): The imputation method.

    Returns:
        AbundanceTable: The preprocessed abundance table.
    """
    steps = [
        FilterGlycans(max_na=filter_max_na),
        Impute(method=impute_method),
        Normalize(),
    ]
    pipeline = ProcessingPipeline(steps=steps)
    return pipeline(data)


class ProcessingStep(Protocol):
    """The protocol for processing steps."""

    def __call__(self, data: AbundanceTable) -> AbundanceTable: ...


@define
class ProcessingPipeline:
    """The pipeline for processing the abundance table."""

    _steps: list[ProcessingStep]

    def __call__(self, data: AbundanceTable) -> AbundanceTable:
        for step in self._steps:
            data = step(data)
        return data


@define
class FilterGlycans(ProcessingStep):
    """Filter glycans with too many missing values.

    Args:
        max_na (float): The maximum proportion of missing values allowed for a glycan.
    """

    max_na: float

    def __call__(self, data: AbundanceTable) -> AbundanceTable:
        to_keep_mask = data.isna().mean() <= self.max_na
        return AbundanceTable(data.loc[:, to_keep_mask])


@define
class Impute(ProcessingStep):
    """Impute the missing values of the abundance table.

    The following imputation methods supported:
        - "zero": Replace the missing values with 0.
        - "min": Replace the missing values with the minimum value of the corresponding glycan.
        - "lod": Replace the missing values with 1/5 of the minimum value of the corresponding
            glycan.
        - "mean": Replace the missing values with the mean value of the corresponding glycan.
        - "median": Replace the missing values with the median value of the corresponding glycan.

    Args:
        method (str): The imputation method.
            Can be "zero", "min", "lod", "mean", "median".
    """

    method: Literal["zero", "min", "lod", "mean", "median"] = field(
        validator=in_(["zero", "min", "lod", "mean", "median"])
    )

    def __call__(self, data: AbundanceTable) -> AbundanceTable:
        if self.method == "zero":
            imputed_df = data.fillna(0)
        elif self.method == "min":
            imputed_df = data.fillna(data.min())
        elif self.method == "lod":
            imputed_df = data.fillna(data.min() / 5)
        elif self.method == "mean":
            imputed_df = data.fillna(data.mean())
        elif self.method == "median":
            imputed_df = data.fillna(data.median())
        return AbundanceTable(imputed_df)


@define
class Normalize(ProcessingStep):
    """Normalize the abundance table by dividing the sum of each sample."""

    def __call__(self, data: AbundanceTable) -> AbundanceTable:
        row_sums = data.sum(axis=1)
        normalized_df = data.div(row_sums, axis=0)
        return AbundanceTable(normalized_df)
