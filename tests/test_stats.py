import pandas as pd
import pytest

from glytrait import stats
from glytrait.exception import HypothesisTestingError


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
        "p-val-adjusted",
    ]
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()


def test_anova():
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 2, 3, 4, 5],
            "trait_2": [6, 7, 8, 9, 10],
            "trait_3": [11, 12, 13, 14, 15],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_2", "group_2", "group_3"])
    result = stats.anova(trait_df, groups)
    expected_columns = [
        "ddof1",
        "ddof2",
        "F",
        "p-val",
        "np2",
        "reject",
        "p-val-adjusted",
        "posthoc",
    ]
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()


def test_auto_hypothesis_test_2_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4])
    ttest_mock = mocker.patch("glytrait.stats.ttest")
    groups = pd.Series(["a", "a", "b", "b"], index=[1, 2, 3, 4])
    stats.auto_hypothesis_test(trait_df, groups)
    ttest_mock.assert_called_once_with(trait_df, groups)


def test_auto_hypothesis_test_3_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4, 5, 6])
    anova_mock = mocker.patch("glytrait.stats.anova")
    groups = pd.Series(["a", "a", "b", "b", "c", "c"], index=[1, 2, 3, 4, 5, 6])
    stats.auto_hypothesis_test(trait_df, groups)
    anova_mock.assert_called_once_with(trait_df, groups)


def test_auto_hypothesis_test_1_group(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3])
    groups = pd.Series(["a", "a", "a"], index=[1, 2, 3])
    with pytest.raises(HypothesisTestingError) as excinfo:
        stats.auto_hypothesis_test(trait_df, groups)
    assert "Only one group is provided." in str(excinfo.value)


def test_auto_hypothesis_test_index_not_same(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3])
    groups = pd.Series(["a", "a", "b"], index=[1, 2, 4])
    with pytest.raises(HypothesisTestingError) as excinfo:
        stats.auto_hypothesis_test(trait_df, groups)
    assert "The index of groups and trait_df must be the same." in str(excinfo.value)


def test_auto_hypothesis_test_index_not_same_order(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 3, 2])
    ttest_mock = mocker.patch("glytrait.stats.ttest")
    groups = pd.Series(["a", "a", "b"], index=[1, 2, 3])
    stats.auto_hypothesis_test(trait_df, groups)
    ttest_mock.assert_called_once_with(trait_df, groups)
