import pandas as pd
import pytest

from glytrait import meta_property as mp
from .glycoct import *


# !!! All tests on meta-properties should use this function for conciseness. !!!
def _test_meta_property(meta_property_type, string, expected, make_glycan):
    """Test a meta-property on a glycan string.

    Args:
        meta_property_type (Type[mp.MetaProperty]): Meta property type to test.
        string (str): Glycan string. Either GlycoCT or composition.
        expected (Any): Expected result.
        make_glycan (Callable[[str], glyc.Structure]): Factory fixture to make a
            glycan from a string. Either `make_structure` or `make_composition`.
            See `conftest.py` for these two fixtures.
    """
    mp_ = meta_property_type()
    glycan = make_glycan(string)
    assert mp_(glycan) == expected


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_2, 2),
        (test_glycoct_3, 0),  # High-mannose
        (test_glycoct_4, 0),  # Hybrid
        (test_glycoct_6, 1),
        (test_glycoct_8, 3),
        (test_glycoct_14, 4),
    ]
)
def test_count_antenna_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.CountAntennaMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 3),
    ],
)
def test_count_fuc_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountFucMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 1),
        ("H5N4", 0),
        ("H5N4F2E1", 2),
    ],
)
def test_count_fuc_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountFucMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 1),
    ],
)
def test_count_core_fuc_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.CountCoreFucMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 0),
        (test_glycoct_9, 2),
    ],
)
def test_count_antennary_fuc_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.CountAntennaryFucMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_2, 1),
        (test_glycoct_3, 0),
        (test_glycoct_4, 1),
        (test_glycoct_10, 2),
    ],
)
def test_count_sia_struc_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.CountSiaMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H5N4", 0),
        ("H5N4F2E1", 1),
    ],
)
def test_count_sia_comp_mp(comp, expected, make_composition):
    _test_meta_property(mp.CountSiaMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, "complex"),
        (test_glycoct_2, "complex"),
        (test_glycoct_3, "high_mannose"),
        (test_glycoct_4, "hybrid"),
        (test_glycoct_5, "complex"),
        (test_glycoct_6, "complex"),
        (test_glycoct_11, "complex"),
    ]
)
def test_glycan_type_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.GlycanTypeMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 1),
        (test_glycoct_2, 1),
        (test_glycoct_3, 1),
    ],
)
def test_wildcard(glycoct, expected, make_structure):
    _test_meta_property(mp.Wildcard, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, True),  # Complex
        (test_glycoct_3, False),  # High-mannose
        (
                test_glycoct_4,
                False,
        ),  # Hybrid
    ],
)
def test_is_complex(glycoct, expected, make_structure):
    _test_meta_property(mp.IsComplex, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),  # Complex
        (test_glycoct_3, True),  # High-mannose
        (
                test_glycoct_4,
                False,
        ),  # Hybrid
    ],
)
def test_is_high_mannose(glycoct, expected, make_structure):
    _test_meta_property(mp.IsHighMannose, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),  # Complex
        (test_glycoct_3, False),  # High-mannose
        (
                test_glycoct_4,
                True,
        ),  # Hybrid
    ],
)
def test_is_hybrid(glycoct, expected, make_structure):
    _test_meta_property(mp.IsHybrid, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_2, True),
        (test_glycoct_3, False),
        (test_glycoct_4, False),
    ],
)
def test_is_bisecting(glycoct, expected, make_structure):
    _test_meta_property(mp.IsBisecting, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_2, False),
        (test_glycoct_3, False),
        (test_glycoct_4, False),
        (test_glycoct_5, False),
        (test_glycoct_6, True),
    ],
)
def test_is_1_antennary(glycoct, expected, make_structure):
    _test_meta_property(mp.Is1Antennary, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, True),
        (test_glycoct_2, True),
        (test_glycoct_3, False),
        (test_glycoct_4, False),
        (test_glycoct_5, False),
        (test_glycoct_6, False),
    ],
)
def test_is_2_antennary(glycoct, expected, make_structure):
    _test_meta_property(mp.Is2Antennary, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_8, True),
        (test_glycoct_14, False),
    ],
)
def test_is_3_antennary(glycoct, expected, make_structure):
    _test_meta_property(mp.Is3Antennary, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_8, False),
        (test_glycoct_14, True),
    ],
)
def test_is_4_antennary(glycoct, expected, make_structure):
    _test_meta_property(mp.Is4Antennary, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_6, 1),
        (test_glycoct_1, 2),
        (test_glycoct_8, 3),
        (test_glycoct_14, 4),
        (test_glycoct_3, 0),
        (test_glycoct_4, 0),
    ],
)
def test_total_antenna(glycoct, expected, make_structure):
    _test_meta_property(mp.TotalAntenna, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 1),
    ],
)
def test_core_fuc(glycoct, expected, make_structure):
    _test_meta_property(mp.CoreFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 0),
        (test_glycoct_9, 2),
    ],
)
def test_antennary_fuc(glycoct, expected, make_structure):
    _test_meta_property(mp.AntennaryFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 3),
    ],
)
def test_total_fuc_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.TotalFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 1),
        ("H5N4", 0),
        ("H5N4F2E1", 2),
    ],
)
def test_total_fuc_comp(comp, expected, make_composition):
    _test_meta_property(mp.TotalFuc, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_2, True),
        (test_glycoct_9, True),
    ],
)
def test_has_fuc_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.HasFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", True),
        ("H5N4", False),
        ("H5N4F2E1", True),
    ],
)
def test_has_fuc_comp(comp, expected, make_composition):
    _test_meta_property(mp.HasFuc, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, True),
        (test_glycoct_2, False),
        (test_glycoct_9, False),
    ],
)
def test_no_fuc_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.NoFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, True),
        (test_glycoct_2, False),
        (test_glycoct_9, False),
    ],
)
def test_no_fuc_comp(glycoct, expected, make_structure):
    _test_meta_property(mp.NoFuc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_2, 1),
        (test_glycoct_3, 0),
        (test_glycoct_4, 1),
        (test_glycoct_10, 2),
    ],
)
def test_total_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.TotalSia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H5N4", 0),
        ("H5N4F2E1", 1),
    ],
)
def test_total_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.TotalSia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, True),
        (test_glycoct_2, True),
        (test_glycoct_3, False),
        (test_glycoct_4, True),
        (test_glycoct_10, True),
    ],
)
def test_has_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.HasSia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", True),
        ("H5N4", False),
        ("H5N4F2E1", True),
    ],
)
def test_has_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.HasSia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_2, False),
        (test_glycoct_3, True),
        (test_glycoct_4, False),
        (test_glycoct_10, False),
    ],
)
def test_no_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.NoSia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", False),
        ("H5N4", True),
        ("H5N4F2E1", False),
    ],
)
def test_no_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.NoSia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 3),
        (test_glycoct_3, 5),
        (test_glycoct_5, 3),
    ],
)
def test_total_man_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.TotalMan, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 3),
        ("H6N5", 3),
        ("H5N2", 5),
    ],
)
def test_total_man_comp(comp, expected, make_composition):
    _test_meta_property(mp.TotalMan, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_3, 0),
        (test_glycoct_4, 1),
    ],
)
def test_total_gal_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.TotalGal, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H6N5", 3),
        ("H5N2", 0),
    ],
)
def test_total_gal_comp(comp, expected, make_composition):
    _test_meta_property(mp.TotalGal, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_14, True),
    ],
)
def test_has_poly_lacnac(glycoct, expected, make_structure):
    _test_meta_property(mp.HasPolyLacNAc, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 0),
        (test_glycoct_11, 2),
        (test_glycoct_12, 1),
    ],
)
def test_a23_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.A23Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", 0),
        ("H5N4F1L2", 2),
        ("H5N4F1L1E1", 1),
        ("H5N4F1E2", 0),
    ],
)
def test_a23_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.A23Sia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 2),
        (test_glycoct_11, 0),
        (test_glycoct_12, 1),
    ],
)
def test_a26_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.A26Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", 0),
        ("H5N4F1L2", 0),
        ("H5N4F1L1E1", 1),
        ("H5N4F1E2", 2),
    ],
)
def test_a26_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.A26Sia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, False),
        (test_glycoct_11, True),
        (test_glycoct_12, True),
    ],
)
def test_hasa23_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.HasA23Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", False),
        ("H5N4F1L2", True),
        ("H5N4F1L1E1", True),
        ("H5N4F1E2", False),
    ],
)
def test_hasa23_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.HasA23Sia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, True),
        (test_glycoct_11, False),
        (test_glycoct_12, True),
    ],
)
def test_hasa26_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.HasA26Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", False),
        ("H5N4F1L2", False),
        ("H5N4F1L1E1", True),
        ("H5N4F1E2", True),
    ],
)
def test_hasa26_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.HasA26Sia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, True),
        (test_glycoct_11, False),
        (test_glycoct_12, False),
    ],
)
def test_noa23_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.NoA23Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", True),
        ("H5N4F1L2", False),
        ("H5N4F1L1E1", False),
        ("H5N4F1E2", True),
    ],
)
def test_noa23_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.NoA23Sia, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, False),
        (test_glycoct_11, True),
        (test_glycoct_12, False),
    ],
)
def test_noa26_sia_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.NoA26Sia, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", True),
        ("H5N4F1L2", True),
        ("H5N4F1L1E1", False),
        ("H5N4F1E2", False),
    ],
)
def test_noa26_sia_comp(comp, expected, make_composition):
    _test_meta_property(mp.NoA26Sia, comp, expected, make_composition)


@pytest.mark.parametrize("comp, expected", [("H5N4F1", False), ("H6N5F1S1", True)])
def test_is_high_branching(comp, expected, make_composition):
    _test_meta_property(mp.IsHighBranching, comp, expected, make_composition)


@pytest.mark.parametrize("comp, expected", [("H5N4F1", True), ("H6N5F1S1", False)])
def test_is_low_branching(comp, expected, make_composition):
    _test_meta_property(mp.IsLowBranching, comp, expected, make_composition)


def test_build_meta_property_table(mocker, make_structure):
    available_mp_mock = mocker.patch(
        "glytrait.meta_property.available_meta_properties",
        return_value={".": mp.Wildcard, "isComplex": mp.IsComplex},
        autospec=True,
    )
    glycans = {
        "glycan1": make_structure(test_glycoct_1),
        "glycan2": make_structure(test_glycoct_2),
        "glycan3": make_structure(test_glycoct_3),
    }

    result = mp.build_meta_property_table(glycans, mode="structure", sia_linkage=False)
    expected = pd.DataFrame(
        {
            ".": [1.0, 1.0, 1.0],
            "isComplex": [True, True, False],
        },
        index=list(glycans.keys()),
    )

    pd.testing.assert_frame_equal(result, expected)
    available_mp_mock.assert_called_once()
