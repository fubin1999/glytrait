import numpy as np
import pandas as pd

from glytrait.stat import t_test, anova


def test_t_test():
    trait_df = pd.DataFrame(
        {
            "Trait1": [1, 2, 3, 4, 5, 6],
            "Trait2": [2, 3, 4, np.nan, 6, 7],
            "Trait3": [3, 4, 5, 6, 7, 8],
        },
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    groups = pd.Series(
        ["A", "A", "A", "B", "B", "B"],
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    result = t_test(trait_df, groups)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 3


def test_anova():
    trait_df = pd.DataFrame(
        {
            "Trait1": list(range(1, 10)),
            "Trait2": list(range(2, 11)),
            "Trait3": list(range(3, 12)),
        },
        index=[f"S{i}" for i in range(1, 10)],
    )
    groups = pd.Series(
        ["A", "A", "A", "B", "B", "B", "C", "C", "C"],
        index=[f"S{i}" for i in range(1, 10)],
    )
    anova_result, post_hoc_result = anova(trait_df, groups)
    assert isinstance(anova_result, pd.DataFrame)
    assert isinstance(post_hoc_result, pd.DataFrame)
    assert anova_result.shape[0] == 3
