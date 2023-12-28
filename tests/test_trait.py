import pandas as pd

import glytrait.formula
import glytrait.meta_property
import glytrait.post_filtering
from glytrait import trait


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
        dtype=float,
    )

    meta_prop_df = pd.DataFrame(
        {
            ".": [1, 1, 1],
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
