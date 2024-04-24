import pandas as pd
import pytest

from glytrait.stat import auto_test


class TestAutoTest:

    @pytest.fixture
    def trait_df(self):
        return pd.DataFrame(
            {"Trait1": [1, 2, 3, 4, 5, 6], "Trait2": [2, 3, 4, 5, 6, 7]},
            index=[f"S{i}" for i in range(1, 7)],
            dtype=float,
        )

    def test_two_groups(self, trait_df, mocker):
        mocker.patch("glytrait.stat.t_test", return_value="t_test_result")
        groups = pd.Series(
            ["A", "A", "A", "B", "B", "B"],
            index=[f"S{i}" for i in range(1, 7)],
            name="Group",
        )
        results = auto_test(trait_df, groups)
        assert results == {"t_test": "t_test_result"}

    def test_more_than_two_groups(self, trait_df, mocker):
        mocker.patch(
            "glytrait.stat.anova", return_value=("anova_result", "post_hoc_result")
        )
        groups = pd.Series(
            ["A", "A", "B", "B", "C", "C"],
            index=[f"S{i}" for i in range(1, 7)],
            name="Group",
        )
        results = auto_test(trait_df, groups)
        assert results == {"anova": "anova_result", "post_hoc": "post_hoc_result"}

    def test_one_group(self, trait_df):
        groups = pd.Series(
            ["A", "A", "A", "A", "A", "A"],
            index=[f"S{i}" for i in range(1, 7)],
            name="Group",
        )
        with pytest.raises(ValueError):
            auto_test(trait_df, groups)
