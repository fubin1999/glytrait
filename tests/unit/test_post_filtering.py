import numpy as np
import pandas as pd
import pytest
from attr import define

import glytrait.post_filtering as pf
import glytrait.formula as fml


def test_filter_invalid():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [np.nan, np.nan, np.nan],
            "C": [1, 1, np.nan],
            "D": [0, 0, 0],
            "E": [1, 1, 1],
            "F": [0.5, 0.5, 0.5],
        }
    )
    result_df = pf.filter_invalid(df)
    expected_df = pd.DataFrame(
        {
            "A": [1, 2, 3],
        }
    )
    pd.testing.assert_frame_equal(result_df, expected_df, check_dtype=False)


def test_relationship_matrix(mocker):
    def mock_is_child_of(trait1, trai2):
        if trait1.name == "trait1":
            return False
        elif trait1.name == "trait2":
            return trai2.name == "trait1"
        elif trait1.name == "trait3":
            return trai2.name in ["trait1", "trait2"]

    trait1 = mocker.Mock()
    trait1.name = "trait1"
    trait2 = mocker.Mock()
    trait2.name = "trait2"
    trait3 = mocker.Mock()
    trait3.name = "trait3"

    mocker.patch("glytrait.post_filtering._is_child_of", side_effect=mock_is_child_of)
    result = pf._relationship_matrix(
        ["trait1", "trait2", "trait3"], [trait1, trait2, trait3]
    )
    expected = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
    assert np.array_equal(result, expected)


@pytest.mark.parametrize(
    "threshold, expected",
    [
        (
            0.5,
            np.array(
                [
                    [1, 1, 0, 1],
                    [1, 1, 0, 1],
                    [0, 0, 1, 0],
                    [1, 1, 0, 1],
                ]
            ),
        ),
        (
            -1,
            np.ones((4, 4), dtype=int),
        ),
        (
            1,
            np.array(
                [
                    [1, 1, 0, 0],
                    [1, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1],
                ]
            ),
        ),
    ],
)
def test_correlation_matrix(threshold, expected):
    df = pd.DataFrame(
        {
            "trait1": [1, 2, 3, 4, 5],
            "trait2": [1, 2, 3, 4, 5],  # 完全与trait1相关
            "trait3": [5, 4, 3, 2, 1],  # 完全与trait1负相关
            "trait4": [2, 3, 1, 5, 4],  # 与trait1的相关性为0.6
        }
    )
    result = pf._correlation_matrix(df, threshold, method="pearson")
    assert np.array_equal(result, expected)


def test_filter_colinearity(mocker):
    trait_table = pd.DataFrame(
        {  # This df is just a place-holder. The values are not important.
            "trait1": [1, 2, 3, 4, 5],
            "trait2": [1, 2, 3, 4, 5],
            "trait3": [5, 4, 3, 2, 1],
            "trait4": [2, 3, 1, 5, 4],
        }
    )
    relationship_matrix_mock = mocker.patch(
        "glytrait.post_filtering._relationship_matrix",
        return_value=np.array(
            [
                [0, 1, 0, 1],
                [0, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ]
        ),
        autospec=True,
    )
    correlation_matrix_mock = mocker.patch(
        "glytrait.post_filtering._correlation_matrix",
        return_value=np.array(
            [
                [1, 1, 0, 1],
                [1, 1, 0, 1],
                [0, 0, 1, 0],
                [1, 1, 0, 1],
            ]
        ),
        autospec=True,
    )
    result_df = pf.filter_colinearity("formulas", trait_table, 0.5, "pearson")
    expected_df = pd.DataFrame(
        {
            "trait2": [1, 2, 3, 4, 5],
            "trait3": [5, 4, 3, 2, 1],
            "trait4": [2, 3, 1, 5, 4],
        }
    )

    pd.testing.assert_frame_equal(result_df, expected_df, check_dtype=False)
    relationship_matrix_mock.assert_called_once_with(
        ["trait1", "trait2", "trait3", "trait4"], "formulas"
    )
    correlation_matrix_mock.assert_called_once_with(trait_table, 0.5, "pearson")


@pytest.mark.parametrize(
    "trait1, trait2, expected",
    [
        ("A2G", "CG", True),
        ("A2Fa", "CFa", True),
        ("A2Fc", "CFc", True),
        ("A2S", "CS", True),
        ("A2E", "CE", True),
        ("A2SG", "A2G", True),
        ("A2FSG", "A2FG", True),
        ("A2FSG", "A2SG", True),
        ("A2F0G", "A2FG", False),
        ("A2FSG", "A2G", False),
    ],
)
def test_is_child_of(trait1, trait2, expected):
    formulas = fml.load_formulas("structure", sia_linkage=True)
    formula_map = {f.name: f for f in formulas}
    formulas1 = formula_map[trait1]
    formulas2 = formula_map[trait2]
    assert pf._is_child_of(formulas1, formulas2) == expected
