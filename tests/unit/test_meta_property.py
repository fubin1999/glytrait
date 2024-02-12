import pandas as pd
import pytest

from glytrait import meta_property as mp
from .glycoct import *


def test_register():
    @mp._register_mp
    class SomeMP(mp.MetaProperty):
        name = "some_mp"
        supported_mode = "structure"

        def __call__(self, glycan) -> int:
            return 0

    assert mp._mp_objects["some_mp"] == SomeMP()


def test_register_duplicate():
    @mp._register_mp
    class MP1:
        name = "mp1"

    with pytest.raises(ValueError):
        @mp._register_mp
        class MP2:
            name = "mp1"


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


def test_mp_return_type():
    """Test the `return_type` property of a meta-property class."""
    class SomeMP(mp.MetaProperty):
        name = "some_mp"
        supported_mode = "structure"

        def __call__(self, glycan) -> int:
            return 0

    assert SomeMP().return_type is int


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
        (test_glycoct_1, False),
        (test_glycoct_2, True),
        (test_glycoct_3, False),
        (test_glycoct_4, False),
        (test_glycoct_8, True),
    ]
)
def test_bisection_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.BisectionMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 3),
        (test_glycoct_3, 5),
        (test_glycoct_5, 3),
    ],
)
def test_count_man_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountManMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 3),
        ("H6N5", 3),
        ("H5N2", 5),
    ],
)
def test_count_man_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountManMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_3, 0),
        (test_glycoct_4, 1),
    ],
)
def test_count_gal_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountGalMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H6N5", 3),
        ("H5N2", 0),
    ],
)
def test_count_gal_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountGalMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 4),
        (test_glycoct_2, 5),
        (test_glycoct_3, 2),
    ]
)
def test_count_glcnac_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountGlcNAcMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 4),
        ("H6N5", 5),
        ("H5N2", 2),
    ]
)
def test_count_glcnac_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountGlcNAcMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_14, True),
    ],
)
def test_has_poly_lacnac_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.HasPolyLacNAcMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 0),
        (test_glycoct_11, 2),
        (test_glycoct_12, 1),
    ],
)
def test_count_a23_sia_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountA23SiaMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", 0),
        ("H5N4F1L2", 2),
        ("H5N4F1L1E1", 1),
        ("H5N4F1E2", 0),
    ],
)
def test_count_a23_sia_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountA23SiaMP, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 2),
        (test_glycoct_11, 0),
        (test_glycoct_12, 1),
    ],
)
def test_count_a26_sia_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.CountA26SiaMP, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1", 0),
        ("H5N4F1L2", 0),
        ("H5N4F1L1E1", 1),
        ("H5N4F1E2", 2),
    ],
)
def test_count_a26_sia_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.CountA26SiaMP, comp, expected, make_composition)
