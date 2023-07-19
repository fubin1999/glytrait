import pandas as pd
import pytest

import glytrait.meta_property
from glytrait import glycan as glyc
from tests import glycoct


def test_build_meta_property_table_struc_no_sia_linkage(make_glycan):
    glycoct_strings = {
        k: v for k, v in glycoct.__dict__.items() if k.startswith("test_glycoct")
    }
    del glycoct_strings["test_glycoct_13"]  # an O-glycan

    glycan_ids = list(glycoct_strings.keys())
    glycans = [make_glycan(s) for s in glycoct_strings.values()]

    result = glytrait.meta_property.build_meta_property_table(
        glycan_ids, glycans, sia_linkage=False, mode="structure"
    )
    expected = pd.DataFrame(
        {
            "isComplex": [1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
            "isHighMannose": [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            "isHybrid": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            "isBisecting": [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            "is1Antennary": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            "is2Antennary": [1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0],
            "is3Antennary": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            "is4Antennary": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            "totalAntenna": [2, 2, 0, 0, 0, 1, 0, 3, 2, 2, 2, 2, 4],
            "coreFuc": [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
            "antennaryFuc": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
            "hasAntennaryFuc": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            "totalFuc": [0, 1, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1],
            "hasFuc": [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
            "noFuc": [1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0],
            "totalSia": [2, 1, 0, 1, 0, 1, 0, 1, 1, 2, 2, 2, 4],
            "hasSia": [1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1],
            "noSia": [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
            "totalMan": [3, 3, 5, 5, 3, 3, 6, 3, 3, 3, 3, 3, 3],
            "totalGal": [2, 2, 0, 1, 0, 1, 0, 3, 1, 2, 2, 2, 5],
            "hasPolyLacNAc": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        },
        index=glycan_ids,
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)


def test_build_meta_property_table_struc_sia_linkage(make_glycan):
    glycoct_strings = {
        k: v
        for k, v in glycoct.__dict__.items()
        if k in ("test_glycoct_10", "test_glycoct_11", "test_glycoct_12")
    }

    glycan_ids = list(glycoct_strings.keys())
    glycans = [make_glycan(s) for s in glycoct_strings.values()]

    result = glytrait.meta_property.build_meta_property_table(
        glycan_ids, glycans, sia_linkage=True, mode="structure"
    )
    partial_expected = pd.DataFrame(
        {
            "a23Sia": [0, 2, 1],
            "a26Sia": [2, 0, 1],
            "hasa23Sia": [0, 1, 1],
            "hasa26Sia": [1, 0, 1],
            "noa23Sia": [1, 0, 0],
            "noa26Sia": [0, 1, 0],
        },
        index=glycan_ids,
    )
    assert len(result.columns) == 27
    partial_result = result[partial_expected.columns]
    pd.testing.assert_frame_equal(partial_result, partial_expected, check_dtype=False)


def test_build_meta_property_table_comp_no_sia_linkage():
    comps = ["H5N4", "H7N2", "H5N4F1S1", "H5N5F1S2", "H5N4S1", "H5N4F2"]
    glycans = [glyc.Composition.from_string(s, sia_linkage=False) for s in comps]
    result = glytrait.meta_property.build_meta_property_table(
        comps, glycans, mode="composition", sia_linkage=False
    )
    expected = pd.DataFrame(
        {
            "isHighBranching": [0, 0, 0, 1, 0, 0],
            "isLowBranching": [1, 1, 1, 0, 1, 1],
            "totalSia": [0, 0, 1, 2, 1, 0],
            "totalFuc": [0, 0, 1, 1, 0, 2],
            "totalGal": [2, 0, 2, 2, 2, 2],
            "hasSia": [0, 0, 1, 1, 1, 0],
            "hasFuc": [0, 0, 1, 1, 0, 1],
            "hasGal": [1, 0, 1, 1, 1, 1],
            "noSia": [1, 1, 0, 0, 0, 1],
            "noFuc": [1, 1, 0, 0, 1, 0],
            "noGal": [0, 1, 0, 0, 0, 0],
        },
        index=comps,
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)


def test_build_meta_property_table_comp_sia_linkage():
    comps = ["H5N4", "H5N4E1", "H5N4L1", "H5N4E1L1", "H5N4E1L2"]
    glycans = [glyc.Composition.from_string(s, sia_linkage=True) for s in comps]
    result = glytrait.meta_property.build_meta_property_table(
        comps, glycans, mode="composition", sia_linkage=True
    )
    partial_expected = pd.DataFrame(
        {
            "totalSia": [0, 1, 1, 2, 3],
            "hasSia": [0, 1, 1, 1, 1],
            "noSia": [1, 0, 0, 0, 0],
            "a23Sia": [0, 0, 1, 1, 2],
            "a26Sia": [0, 1, 0, 1, 1],
            "hasa23Sia": [0, 0, 1, 1, 1],
            "hasa26Sia": [0, 1, 0, 1, 1],
            "noa23Sia": [1, 1, 0, 0, 0],
            "noa26Sia": [1, 0, 1, 0, 0],
        },
        index=comps,
    )
    assert len(result.columns) == 17
    partial_result = result[partial_expected.columns]
    pd.testing.assert_frame_equal(partial_result, partial_expected, check_dtype=False)
