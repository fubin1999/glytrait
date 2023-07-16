import warnings

import pandas as pd
import pingouin as pg
from sklearn.metrics import roc_auc_score

from glytrait.exception import HypothesisTestingError
from glytrait.trait import filter_derived_trait


def mwu(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Perform Mann-Whitney U Test for two groups.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.
            Only two groups are allowed.

    Returns:
        pd.DataFrame: Dataframe containing the Mann-Whitney U Test results,
            with trait names as index.
    """
    ttest_results = []
    group_names = groups.unique()
    group_1_idx = groups == group_names[0]
    group_2_idx = groups == group_names[1]
    for col in trait_df.columns:
        ttest_result = pg.mwu(
            trait_df.loc[group_1_idx, col], trait_df.loc[group_2_idx, col]
        )
        ttest_result.index = [col]
        ttest_results.append(ttest_result)
    result_df = pd.concat(ttest_results, axis=0)
    _, result_df["p-val-adjusted"] = pg.multicomp(result_df["p-val"], method="fdr_bh")
    result_df = result_df.drop(["alternative", "RBC", "CLES"], axis=1)
    return result_df


def kruskal(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Perform Kruskal-Wallis H-test for multiple groups.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.

    Returns:
        pd.DataFrame: Dataframe containing the Kruskal-Wallis H-test results,
            with trait names as index.
    """
    anove_results = []
    trait_names = trait_df.columns
    groups = groups.rename("group")
    trait_df = trait_df.merge(groups, left_index=True, right_index=True)
    for trait in trait_names:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            anova_result = pg.kruskal(data=trait_df, dv=trait, between="group")
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
            posthoc_result = pg.pairwise_tests(
                data=trait_df, dv=trait, between="group", parametric=False
            )
            posthoc_result = posthoc_result[posthoc_result["p-unc"] < 0.05]
            result_df.loc[trait, "posthoc"] = ",".join(
                posthoc_result["A"] + "-" + posthoc_result["B"]
            )

    result_df = result_df.drop(["reject", "ddof1"], axis=1)
    return result_df


def auto_hypothesis_test(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Automatically perform hypothesis test for the trait data.

    If `groups` has two unique values, Mann-Whitney U Test will be performed.
    Otherwise, Kruskal-Wallis H-test will be performed.

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
    trait_df = filter_derived_trait(trait_df)
    if set(groups.index) != set(trait_df.index):
        raise HypothesisTestingError(
            "The index of groups and trait_df must be the same."
        )
    if groups.nunique() == 1:
        raise HypothesisTestingError("Only one group is provided.")
    elif groups.nunique() == 2:
        return mwu(trait_df, groups)
    else:
        return kruskal(trait_df, groups)


def calcu_roc_auc(trait_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Perform ROC AUC test for the trait data.

    Args:
        trait_df (pd.DataFrame): Dataframe containing the trait data.
        groups (pd.Series): Series containing the group information, with the same index as df.

    Returns:
        pd.DataFrame: Dataframe containing the ROC AUC test results with two columns:
            "trait" and "roc_auc".
    """
    roc_auc_results = []
    for trait in trait_df.columns:
        roc_auc = roc_auc_score(groups, trait_df[trait])
        roc_auc_results.append((trait, roc_auc))
    result_df = pd.DataFrame.from_records(roc_auc_results, columns=["trait", "ROC AUC"])
    return result_df
