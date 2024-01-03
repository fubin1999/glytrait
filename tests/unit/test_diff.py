import numpy as np
import pandas as pd
import pytest

import glytrait.diff


def test_mwu():
    trait_df = pd.DataFrame(
        {
            "trait_1": [1, 2, 3, 4, 5],
            "trait_2": [6, 7, 8, 9, 10],
            "trait_3": [11, 12, 13, 14, 15],
        }
    )
    groups = pd.Series(["group_1", "group_1", "group_1", "group_2", "group_2"])
    result = glytrait.diff.mwu(trait_df, groups)
    expected_columns = [
        "U-val",
        "p-val",
        "p-val-adjusted",
        "RBC",
        "CLES",
    ]
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()
    assert result.index.name == "trait"


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

    result = glytrait.diff.kruskal(trait_df, groups)
    expected_columns = [
        "H",
        "p-val",
        "p-val-adjusted",
        "posthoc",
    ]
    assert multicomp_mock.call_count == 1
    assert result.columns.tolist() == expected_columns
    assert result.index.tolist() == trait_df.columns.tolist()
    assert result.index.name == "trait"


def test_differential_analysis_2_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4])
    mwu_mock = mocker.patch("glytrait.diff.mwu", autospec=True)
    groups = pd.Series(["a", "a", "b", "b"], index=[1, 2, 3, 4])
    glytrait.diff.differential_analysis(trait_df, groups)
    mwu_mock.assert_called_once_with(trait_df, groups)


def test_differential_analysis_3_groups(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 2, 3, 4, 5, 6])
    kruskal_mock = mocker.patch("glytrait.diff.kruskal", autospec=True)
    groups = pd.Series(["a", "a", "b", "b", "c", "c"], index=[1, 2, 3, 4, 5, 6])
    glytrait.diff.differential_analysis(trait_df, groups)
    kruskal_mock.assert_called_once_with(trait_df, groups)


def test_differential_analysis_index_not_same_order(mocker):
    trait_df = mocker.Mock()
    trait_df.index = pd.Index([1, 3, 2])
    mwu_mock = mocker.patch("glytrait.diff.mwu", autospec=True)
    groups = pd.Series(["a", "a", "b"], index=[1, 2, 3])
    glytrait.diff.differential_analysis(trait_df, groups)
    mwu_mock.assert_called_once_with(trait_df, groups)
