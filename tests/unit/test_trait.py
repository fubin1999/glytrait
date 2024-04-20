import pandas as pd
from attrs import define

from glytrait.trait import calcu_derived_trait


@define
class FakeFormula:

    name: str
    value: float

    def initialize(self, *args, **kwargs):
        return

    def calcu_trait(self, abundance_df: pd.DataFrame):
        return pd.Series(
            [self.value] * len(abundance_df.index),
            index=abundance_df.index,
            name=self.name,
        )


def test_calcu_derived_trait(mocker):
    abund_df = pd.DataFrame(
        {
            "G1": [1, 2, 3],
            "G2": [4, 5, 6],
            "G3": [7, 8, 9],
        },
        index=["S1", "S2", "S3"],
        dtype=float,
    )

    formulas = [FakeFormula("all1", 1), FakeFormula("all2", 2)]

    result = calcu_derived_trait(abund_df, mocker.Mock(), formulas)
    expected = pd.DataFrame(
        {
            "all1": [1, 1, 1],
            "all2": [2, 2, 2],
        },
        index=["S1", "S2", "S3"],
    )
    pd.testing.assert_frame_equal(result, expected)
