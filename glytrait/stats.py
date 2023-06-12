import warnings

import pandas as pd
import pingouin as pg

from .exception import HypothesisTestingError


def ttest(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Perform t-test for two groups.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.
            Only two groups are allowed.

    Returns:
        pd.DataFrame: Dataframe containing the t-test results, with trait names as index.
    """
    ttest_results = []
    group_names = groups.unique()
    group_1_idx = groups == group_names[0]
    group_2_idx = groups == group_names[1]
    for col in trait_df.columns:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            ttest_result = pg.ttest(
                trait_df.loc[group_1_idx, col], trait_df.loc[group_2_idx, col]
            )
        ttest_result.index = [col]
        ttest_results.append(ttest_result)
    result_df = pd.concat(ttest_results, axis=0)
    _, result_df["p-val-adjusted"] = pg.multicomp(result_df["p-val"], method="fdr_bh")
    result_df["CI95%"] = result_df["CI95%"].astype(str)
    return result_df


def anova(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Perform ANOVA for multiple groups.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.

    Returns:
        pd.DataFrame: Dataframe containing the ANOVA results, with trait names as index.
    """
    anove_results = []
    trait_names = trait_df.columns
    groups = groups.rename("group")
    trait_df = trait_df.merge(groups, left_index=True, right_index=True)
    for trait in trait_names:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            anova_result = pg.anova(data=trait_df, dv=trait, between="group")
        anova_result.index = [trait]
        anove_results.append(anova_result)
    result_df = pd.concat(anove_results, axis=0)
    result_df = result_df.drop("Source", axis=1)
    result_df = result_df.rename(columns={"p-unc": "p-val"})
    result_df["reject"], result_df["p-val-adjusted"] = pg.multicomp(
        result_df["p-val"], method="fdr_bh"
    )

    result_df["posthoc"] = ""
    for trait in result_df.index:
        if result_df.loc[trait, "reject"]:
            posthoc_result = pg.pairwise_tukey(data=trait_df, dv=trait, between="group")
            posthoc_result = posthoc_result[posthoc_result["p-tukey"] < 0.05]
            result_df.loc[trait, "posthoc"] = ",".join(
                posthoc_result["A"] + "-" + posthoc_result["B"]
            )

    return result_df


def auto_hypothesis_test(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Automatically perform hypothesis test for the trait data.

    If `groups` has two unique values, t-test will be performed. Otherwise, ANOVA
    will be performed.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.

    Returns:
        pd.DataFrame: Dataframe containing the hypothesis test results,
            with trait names as index.

    Raises:
        HypothesisTestingError: If only one group is provided, or when the index of
            `groups` and `trait_df` are not the same.
    """
    if set(groups.index) != set(trait_df.index):
        raise HypothesisTestingError("The index of groups and trait_df must be the same.")
    if groups.nunique() == 1:
        raise HypothesisTestingError("Only one group is provided.")
    elif groups.nunique() == 2:
        return ttest(trait_df, groups)
    else:
        return anova(trait_df, groups)
