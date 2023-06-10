import warnings

import pandas as pd
import pingouin as pg


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
    return result_df
