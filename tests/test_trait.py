import numpy as np
import pandas as pd
import pytest

import glytrait.formula
import glytrait.meta_property
from glytrait import trait
from glytrait.trait import _relationship_matrix, _correlation_matrix, filter_colinearity


def test_calcu_trait():
    trait_formulas = [
        glytrait.formula.TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["."],
        ),
        glytrait.formula.TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        ),
    ]

    abundance_df = pd.DataFrame(
        {
            "G1": [1, 2, 3],
            "G2": [4, 5, 6],
            "G3": [7, 8, 9],
        },
        index=["S1", "S2", "S3"],
    )

    meta_prop_df = pd.DataFrame(
        {
            "isHighMannose": [True, False, False],
            "isHybrid": [False, False, True],
            "isComplex": [False, True, False],
        },
        index=["G1", "G2", "G3"],
    )

    result = trait.calcu_derived_trait(abundance_df, meta_prop_df, trait_formulas)
    expected = pd.DataFrame(
        {
            "TM": [1 / 12, 2 / 15, 3 / 18],
            "MHy": [1 / 7, 2 / 8, 3 / 9],
        },
        index=["S1", "S2", "S3"],
    )
    pd.testing.assert_frame_equal(result, expected)


def test_filter_derived_trait():
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
    result = trait.filter_invalid(df)
    expected = pd.DataFrame(
        {
            "A": [1, 2, 3],
        }
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)


def test_relationship_matrix(mocker):
    trait1 = mocker.Mock()
    trait1.name = "trait1"
    trait1.is_child_of.return_value = False
    trait2 = mocker.Mock()
    trait2.name = "trait2"
    trait2.is_child_of.side_effect = lambda x: x.name == "trait1"
    trait3 = mocker.Mock()
    trait3.name = "trait3"
    trait3.is_child_of.side_effect = lambda x: x.name in ["trait1", "trait2"]

    result = _relationship_matrix(
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
            np.array(
                [
                    [1, 1, 0, 1],
                    [1, 1, 0, 1],
                    [0, 0, 1, 1],
                    [1, 1, 1, 1],
                ]
            ),
        ),
        (
            1,
            np.zeros((4, 4), dtype=int),
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
    result = _correlation_matrix(df, threshold)
    assert np.array_equal(result, expected)


def test_filter_colinearity(mocker):
    trait_table_mock = mocker.Mock()
    trait_table_mock.columns = ["trait1", "trait2", "trait3", "trait4"]
    relationship_matrix_mock = mocker.patch(
        "glytrait.trait._relationship_matrix",
        return_value=np.array(
            [
                [0, 1, 0, 1],
                [0, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ]
        ),
    )
    correlation_matrix_mock = mocker.patch(
        "glytrait.trait._correlation_matrix",
        return_value=np.array(
            [
                [1, 1, 0, 1],
                [1, 1, 0, 1],
                [0, 0, 1, 0],
                [1, 1, 0, 1],
            ]
        ),
    )
    result = filter_colinearity("formulas", trait_table_mock, 0.5)
    expected = np.array([False, True, True, True])
    assert np.array_equal(result, expected)
    relationship_matrix_mock.assert_called_once_with(
        ["trait1", "trait2", "trait3", "trait4"], "formulas"
    )
    correlation_matrix_mock.assert_called_once_with(trait_table_mock, 0.5)
