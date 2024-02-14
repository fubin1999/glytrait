from itertools import permutations

import pandas as pd
import pytest

from glytrait.formula import TraitFormula
from glytrait.trait import calcu_derived_trait

pytestmark = pytest.mark.skip("`TraitFormula` to be updated.")


@pytest.fixture
def trait_formulas() -> list[TraitFormula]:
    return [
        TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["."],
        ),
        TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        ),
    ]


@pytest.fixture
def abundance_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "G1": [1, 2, 3],
            "G2": [4, 5, 6],
            "G3": [7, 8, 9],
        },
        index=["S1", "S2", "S3"],
        dtype=float,
    )


@pytest.fixture
def meta_property_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            ".": [1, 1, 1],
            "isHighMannose": [True, False, False],
            "isHybrid": [False, False, True],
            "isComplex": [False, True, False],
        },
        index=["G1", "G2", "G3"],
    )


@pytest.fixture
def expected_derived_trait_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "TM": [1 / 12, 2 / 15, 3 / 18],
            "MHy": [1 / 7, 2 / 8, 3 / 9],
        },
        index=["S1", "S2", "S3"],
    )


def test_calcu_trait(
    trait_formulas, abundance_table, meta_property_table, expected_derived_trait_table
):
    result = calcu_derived_trait(abundance_table, meta_property_table, trait_formulas)
    pd.testing.assert_frame_equal(result, expected_derived_trait_table)


def test_calcu_trait_glycan_order_not_the_same(
    trait_formulas, abundance_table, meta_property_table, expected_derived_trait_table
):
    glycans = abundance_table.columns.tolist()
    for perm in permutations(glycans):
        new_mp_table = meta_property_table.reindex(perm)
        result = calcu_derived_trait(abundance_table, new_mp_table, trait_formulas)
        pd.testing.assert_frame_equal(result, expected_derived_trait_table)
