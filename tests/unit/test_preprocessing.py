import numpy as np
import pandas as pd
import pytest
from hypothesis import given, strategies as st
from hypothesis.extra.pandas import columns, data_frames

from glytrait import preprocessing as pp


@given(
    data=data_frames(columns([f"G{i}" for i in range(1, 11)], dtype=float)),
    max_na=st.floats(min_value=0, max_value=1, allow_nan=False),
)
def test_filter(data, max_na):
    data = data.abs()
    filter = pp.FilterGlycans(max_na=max_na)
    result = filter(data)
    # Assert the NA proportions of remaining columns are all less than `max_na`.
    assert np.sum(result.isna().mean() > max_na) == 0


@pytest.mark.parametrize(
    "method, expected",
    [
        (
            "zero",
            pd.DataFrame(
                {
                    "Glycan1": [0.1, 0.2, 0.0, 0.4],
                    "Glycan2": [0.2, 0.0, 0.0, 0.5],
                    "Glycan3": [0.3, 0.4, 0.5, 0.6],
                },
                index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
            ),
        ),
        (
            "min",
            pd.DataFrame(
                {
                    "Glycan1": [0.1, 0.2, 0.1, 0.4],
                    "Glycan2": [0.2, 0.2, 0.2, 0.5],
                    "Glycan3": [0.3, 0.4, 0.5, 0.6],
                },
                index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
            ),
        ),
        (
            "lod",
            pd.DataFrame(
                {
                    "Glycan1": [0.1, 0.2, 0.02, 0.4],
                    "Glycan2": [0.2, 0.04, 0.04, 0.5],
                    "Glycan3": [0.3, 0.4, 0.5, 0.6],
                },
                index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
            ),
        ),
        (
            "mean",
            pd.DataFrame(
                {
                    "Glycan1": [0.1, 0.2, 0.7 / 3, 0.4],
                    "Glycan2": [0.2, 0.35, 0.35, 0.5],
                    "Glycan3": [0.3, 0.4, 0.5, 0.6],
                },
                index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
            ),
        ),
        (
            "median",
            pd.DataFrame(
                {
                    "Glycan1": [0.1, 0.2, 0.2, 0.4],
                    "Glycan2": [0.2, 0.35, 0.35, 0.5],
                    "Glycan3": [0.3, 0.4, 0.5, 0.6],
                },
                index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
            ),
        ),
    ],
)
def test_impute(method, expected):
    data = pd.DataFrame(
        {
            "Glycan1": [0.1, 0.2, None, 0.4],
            "Glycan2": [0.2, None, None, 0.5],
            "Glycan3": [0.3, 0.4, 0.5, 0.6],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="Sample"),
    )
    imputer = pp.Impute(method=method)
    result = imputer(data)
    pd.testing.assert_frame_equal(result, expected)


def test_unknown_impute_method():
    with pytest.raises(ValueError):
        pp.Impute(method="unknown")


@given(
    st.lists(
        st.floats(min_value=0.1, max_value=100, allow_nan=False), min_size=3, max_size=3
    ),
    st.lists(
        st.floats(min_value=0.1, max_value=100, allow_nan=False), min_size=3, max_size=3
    ),
    st.lists(
        st.floats(min_value=0.1, max_value=100, allow_nan=False), min_size=3, max_size=3
    ),
)
def test_normalize(a1, a2, a3):
    data = pd.DataFrame({"G1": a1, "G2": a2, "G3": a3})
    normalizer = pp.Normalize()
    result = normalizer(data)
    np.testing.assert_allclose(result.sum(axis=1).values, np.ones(len(result.index)))
