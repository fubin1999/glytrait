import pandas as pd
from sklearn.metrics import roc_auc_score


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
