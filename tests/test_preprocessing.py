import pandas as pd
import pytest

from glytrait import preprocessing as pp


@pytest.mark.parametrize(
    "method, expected",
    [
        ("zero", 0),
        ("min", 0.1),
        ("lod", 0.02),
        ("mean", 0.2),
        ("median", 0.2),
    ],
)
def test_impute(method, expected):
    df = pd.DataFrame(
        {
            "Glycan1": [0.1, None, 0.3],
            "Glycan2": [0.2, 0.3, 0.4],
            "Glycan3": [0.3, 0.4, 0.5],
        }
    )
    result = pp.impute(df, method)
    assert result.iloc[1, 0] == expected


@pytest.mark.parametrize(
    "ratio, expected",
    [
        (0, ["Glycan2"]),
        (0.3, ["Glycan1", "Glycan2"]),
        (0.5, ["Glycan1", "Glycan2", "Glycan3"]),
        (0.7, ["Glycan1", "Glycan2", "Glycan3", "Glycan4"]),
        (1, ["Glycan1", "Glycan2", "Glycan3", "Glycan4", "Glycan5", "Glycan6"]),
    ],
)
def test_filter_glycans(ratio, expected):
    df = pd.DataFrame(
        {
            "Glycan1": [0.1, None, 0.3, 0.1, 0.2],
            "Glycan2": [0.2, 0.3, 0.4, 0.2, 0.3],
            "Glycan3": [0.3, 0.4, None, None, 0.4],
            "Glycan4": [None, None, None, 0.3, 0.5],
            "Glycan5": [None, None, None, 0.4, None],
            "Glycan6": [None, None, None, None, None],
        }
    )
    result = pp.filter_glycans(df, ratio)
    assert sorted(result.columns.tolist()) == sorted(expected)


def test_normalization():
    df = pd.DataFrame(
        {
            "Glycan1": [0.1, 0.2, 0.3],
            "Glycan2": [0.2, 0.3, 0.4],
            "Glycan3": [0.3, 0.4, 0.5],
        }
    )
    result = pp.normalization(df)
    expected = pd.DataFrame(
        {
            "Glycan1": [0.1 / 0.6, 0.2 / 0.9, 0.3 / 1.2],
            "Glycan2": [0.2 / 0.6, 0.3 / 0.9, 0.4 / 1.2],
            "Glycan3": [0.3 / 0.6, 0.4 / 0.9, 0.5 / 1.2],
        }
    )
    pd.testing.assert_frame_equal(result, expected)
