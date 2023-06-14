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
    ]
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
