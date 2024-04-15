import pandas as pd
import pingouin as pg

from glytrait.data_type import DerivedTraitTable, GroupSeries

__all__ = ["differential_analysis"]


def differential_analysis(
    trait_df: DerivedTraitTable, groups: GroupSeries
) -> pd.DataFrame:
    """Perform hypothesis test for the trait data.

    If `groups` has two unique values, t-test will be performed.
    Otherwise, ANOVA will be performed.

    Args:
        trait_df (DerivedTraitTable): Dataframe containing the trait data.
        groups (GroupSeries): Series containing the group information.

    Returns:
        pd.DataFrame: Dataframe containing the differential analysis results,
            with trait names as index.
    """
    return pd.DataFrame()
