import pandas as pd
import pytest

from glytrait import preprocessing as pp
from glytrait.load_data import GlyTraitInputData


class TestImpute:
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
    def test_impute_logc(self, method, expected):
        df = pd.DataFrame(
            {
                "Glycan1": [0.1, None, 0.3],
                "Glycan2": [0.2, 0.3, 0.4],
                "Glycan3": [0.3, 0.4, 0.5],
            }
        )
        result = pp._impute(df, method)
        assert result.iloc[1, 0] == expected

    def test_impute_step(self):
        df = pd.DataFrame(
            {
                "G1": [1, 2, None],
                "G2": [1, 2, 3],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        data = GlyTraitInputData(
            abundance_table=df, glycans={"G1": "Glycan1", "G2": "Glycan2"}
        )
        step = pp.Impute(method="zero")
        step(data)
        assert data.abundance_table.isna().sum().sum() == 0


class TestFilterGlycans:
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
    def test_filter_glycans_logic(self, ratio, expected):
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
        result = pp._filter_glycans(df, ratio)
        assert sorted(result) == sorted(expected)

    def test_filter_glycans_step(self):
        df = pd.DataFrame(
            {
                "G1": [1, 2, None],
                "G2": [1, 2, 3],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        glycans = {"G1": "Glycan1", "G2": "Glycan2"}
        data = GlyTraitInputData(abundance_table=df, glycans=glycans)
        step = pp.FilterGlycans(max_na=0.2)
        step(data)
        assert data.abundance_table.columns.tolist() == ["G2"]
        assert data.glycans == {"G2": "Glycan2"}


class TestNormalize:
    def test_normalization_logic(self):
        df = pd.DataFrame(
            {
                "Glycan1": [0.1, 0.2, 0.3],
                "Glycan2": [0.2, 0.3, 0.4],
                "Glycan3": [0.3, 0.4, 0.5],
            }
        )
        result = pp._normalization(df)
        expected = pd.DataFrame(
            {
                "Glycan1": [0.1 / 0.6, 0.2 / 0.9, 0.3 / 1.2],
                "Glycan2": [0.2 / 0.6, 0.3 / 0.9, 0.4 / 1.2],
                "Glycan3": [0.3 / 0.6, 0.4 / 0.9, 0.5 / 1.2],
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_normalization_step(self):
        df = pd.DataFrame(
            {
                "G1": [1, 2, 1],
                "G2": [1, 2, 1],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        data = GlyTraitInputData(
            abundance_table=df, glycans={"G1": "Glycan1", "G2": "Glycan2"}
        )
        step = pp.Normalize()
        step(data)
        assert data.abundance_table.sum().sum() == len(df.index)
