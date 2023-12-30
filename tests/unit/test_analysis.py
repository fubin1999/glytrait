import numpy as np
import pandas as pd
import pytest

import glytrait.analysis.diff
import glytrait.analysis.roc
from glytrait.exception import HypothesisTestingError


def test_mwu():
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 2, 3, 4, 5],
            "trait_2": [6, 7, 8, 9, 10],
            "trait_3": [11, 12, 13, 14, 15],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_1", "group_2", "group_2"])
    result = glytrait.analysis.diff.mwu(trait_df, groups)
    expected_columns = [
        "U-val",
        "p-val",
        "p-val-adjusted",
        "CLES",
    ]
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()


def test_kruskal(mocker):
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 1, 3, 4, 5],
            "trait_2": [1, 2, 8, 9, 10],
            "trait_3": [2, 1, 13, 14, 15],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_2", "group_2", "group_2"])

    reject = np.array([True, False, True])
    p_val_adjusted = np.array([0.01, 0.1, 0.01])
    multicomp_mock = mocker.patch(
        "pingouin.multicomp",
        return_value=(reject, p_val_adjusted),
        autospec=True,
    )

    result = glytrait.analysis.diff.kruskal(trait_df, groups)
    expected_columns = [
        "H",
        "p-val",
        "p-val-adjusted",
        "posthoc",
    ]
    assert multicomp_mock.call_count == 1
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()


def test_differential_analysis_2_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4])
    mwu_mock = mocker.patch("glytrait.analysis.diff.mwu", autospec=True)
    groups = pd.Series(["a", "a", "b", "b"], index=[1, 2, 3, 4])
    glytrait.analysis.diff.differential_analysis(trait_df, groups)
    mwu_mock.assert_called_once_with(trait_df, groups)


def test_differential_analysis_3_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4, 5, 6])
    kruskal_mock = mocker.patch("glytrait.analysis.diff.kruskal", autospec=True)
    groups = pd.Series(["a", "a", "b", "b", "c", "c"], index=[1, 2, 3, 4, 5, 6])
    glytrait.analysis.diff.differential_analysis(trait_df, groups)
    kruskal_mock.assert_called_once_with(trait_df, groups)


def test_differential_analysis_1_group(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3])
    groups = pd.Series(["a", "a", "a"], index=[1, 2, 3])
    with pytest.raises(HypothesisTestingError) as excinfo:
        glytrait.analysis.diff.differential_analysis(trait_df, groups)
    assert "Only one group is provided." in str(excinfo.value)


def test_differential_analysis_index_not_same(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3])
    groups = pd.Series(["a", "a", "b"], index=[1, 2, 4])
    with pytest.raises(HypothesisTestingError) as excinfo:
        glytrait.analysis.diff.differential_analysis(trait_df, groups)
    assert "The index of groups and trait_df must be the same." in str(excinfo.value)


def test_differential_analysis_index_not_same_order(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 3, 2])
    mwu_mock = mocker.patch("glytrait.analysis.diff.mwu", autospec=True)
    groups = pd.Series(["a", "a", "b"], index=[1, 2, 3])
    glytrait.analysis.diff.differential_analysis(trait_df, groups)
    mwu_mock.assert_called_once_with(trait_df, groups)


def test_calcu_roc_auc(mocker):
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 2, 3, 4, 5],
            "trait_2": [6, 7, 8, 9, 10],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_1", "group_2", "group_2"])
    roc_auc_score_mock = mocker.patch(
        "glytrait.analysis.roc.roc_auc_score", side_effect=[0.7, 0.8], autospec=True
    )
    result = glytrait.analysis.roc.calcu_roc_auc(trait_df, groups)
    expected = pd.DataFrame({"trait": ["trait_1", "trait_2"], "ROC AUC": [0.7, 0.8]})
    pd.testing.assert_frame_equal(result, expected)

    first_args = [args[0][0] for args in roc_auc_score_mock.call_args_list]
    first_args_expected = [groups, groups]
    for arg1, arg2 in zip(first_args, first_args_expected):
        pd.testing.assert_series_equal(arg1, arg2)

    second_args = [args[0][1] for args in roc_auc_score_mock.call_args_list]
    second_args_expected = [trait_df["trait_1"], trait_df["trait_2"]]
    for arg1, arg2 in zip(second_args, second_args_expected):
        pd.testing.assert_series_equal(arg1, arg2)