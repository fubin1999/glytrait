import pandas as pd

from glytrait import stats


def test_ttest():
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 2, 3, 4, 5],
            "trait_2": [6, 7, 8, 9, 10],
            "trait_3": [11, 12, 13, 14, 15],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_1", "group_2", "group_2"])
    result = stats.ttest(trait_df, groups)
    expected_columns = [
        "T",
        "dof",
        "alternative",
        "p-val",
        "CI95%",
        "cohen-d",
        "BF10",
        "power",
    ]
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()
