import warnings

import pandas as pd
import pingouin as pg

from glytrait.data_type import DerivedTraitTable, GroupSeries

__all__ = ["t_test", "anova"]


def t_test(trait_df: DerivedTraitTable, groups: GroupSeries) -> pd.DataFrame:
    """Perform t-test for the trait data.

    Args:
        trait_df (DerivedTraitTable): Dataframe containing the trait data.
        groups (GroupSeries): Series containing the group information.

    Returns:
        pd.DataFrame: Dataframe containing the t-test results,
            with trait names as index.
    """
    group_names = groups.unique()
    assert len(group_names) == 2, "t-test requires two groups."

    group_data = trait_df.groupby(groups)
    group_1_data = group_data.get_group(group_names[0])
    group_2_data = group_data.get_group(group_names[1])

    results: list[pd.DataFrame] = []
    for trait_name in trait_df.columns:
        t_test_result = pg.ttest(group_1_data[trait_name], group_2_data[trait_name])
        t_test_result.index = [trait_name]
        results.append(t_test_result)
    result_df = pd.concat(results)
    result_df["p-val-adj"] = pg.multicomp(result_df["p-val"].values, method="fdr_bh")[1]
    return result_df


def anova(
    trait_df: DerivedTraitTable, groups: GroupSeries
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform ANOVA for the trait data.

    Args:
        trait_df (DerivedTraitTable): Dataframe containing the trait data.
        groups (GroupSeries): Series containing the group information.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]:
            Tuple containing ANOVA results and post-hoc results.
    """
    trait_names = trait_df.columns.tolist()
    prepared_df = _prepare_for_anova(trait_df, groups)
    anova_result_df = _anova(prepared_df, trait_names)
    sig_traits = anova_result_df[anova_result_df["reject"]].index.tolist()
    post_hoc_result_df = _post_hoc(prepared_df, sig_traits)
    return anova_result_df, post_hoc_result_df


def _prepare_for_anova(
    trait_df: DerivedTraitTable, groups: GroupSeries
) -> pd.DataFrame:
    """Combine the trait data and group information."""
    groups = groups.rename("group")
    trait_df = trait_df.merge(groups, left_index=True, right_index=True)
    return trait_df


def _anova(prepared_df: pd.DataFrame, trait_names: list[str]) -> pd.DataFrame:
    """Perform ANOVA for the trait data.

    Args:
        prepared_df: Dataframe containing the trait data and group information.
        trait_names: List of trait names.
    """
    anove_results: list[pd.DataFrame] = []
    for trait in trait_names:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            anova_result = pg.welch_anova(data=prepared_df, dv=trait, between="group")
        anova_result.index = [trait]
        anove_results.append(anova_result)
    return _tidy_anova_result(pd.concat(anove_results))


def _tidy_anova_result(result_df: pd.DataFrame) -> pd.DataFrame:
    """Tidy the ANOVA result."""
    result_df = result_df.drop("Source", axis=1)
    result_df = result_df.rename(columns={"p-unc": "p-val"})
    result_df["reject"], result_df["p-val-adj"] = pg.multicomp(
        result_df["p-val"], method="fdr_bh"
    )
    return result_df


def _post_hoc(prepared_df: pd.DataFrame, trait_names: list[str]) -> pd.DataFrame:
    """Perform post-hoc test for the trait data."""
    post_hoc_results: list[pd.DataFrame] = []
    for trait_name in trait_names:
        post_hoc_result = pg.pairwise_gameshowell(
            data=prepared_df, dv=trait_name, between="group"
        )
        post_hoc_result.index = [trait_name] * len(post_hoc_result.index)
        post_hoc_results.append(post_hoc_result)
    return pd.concat(post_hoc_results)
