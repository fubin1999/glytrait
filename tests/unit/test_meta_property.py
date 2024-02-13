import pandas as pd
import pytest

from glytrait import meta_property as mp
from .glycoct import *


def test_build_meta_property_table(mocker):
    mocker.patch(
        "glytrait.meta_property.available_meta_properties",
        autospec=True,
        return_value=["mp1", "mp2"],
    )

    mps = {
        "mp1": lambda glycoct: 1,
        "mp2": lambda glycoct: int(glycoct[-1]),
    }
    mocker.patch(
        "glytrait.meta_property.get_meta_property",
        autospec=True,
        side_effect=lambda mp: mps[mp],
    )

    glycans = {"G1": "glycan1", "G2": "glycan2", "G3": "glycan3"}
    result = mp.build_meta_property_table(glycans, "structure", True)
    expected = pd.DataFrame({"mp1": [1, 1, 1], "mp2": [1, 2, 3]}, index=glycans)
    pd.testing.assert_frame_equal(result, expected)


@pytest.fixture
def register_some_mp():
    @mp.mp("some_mp", "both")
    def some_mp(glycan):
        return 0

    yield
    mp._mp_objects.pop("some_mp")


def test_register(register_some_mp):
    assert "some_mp" in mp._mp_objects


def test_register_duplicate(register_some_mp):
    with pytest.raises(ValueError):

        @mp.mp("some_mp", "both")
        def some_mp(glycan):
            return 0


@pytest.mark.parametrize(
    "mode, sia, expected",
    [
        (
            "structure",
            True,
            [
                "type",
                "B",
                "nAnt",
                "nF",
                "nFc",
                "nFa",
                "nS",
                "nM",
                "nG",
                "nN",
                "PL",
                "nL",
                "nE",
            ],
        ),
        (
            "structure",
            False,
            [
                "type",
                "B",
                "nAnt",
                "nF",
                "nFc",
                "nFa",
                "nS",
                "nM",
                "nG",
                "nN",
                "PL",
            ],
        ),
        (
            "composition",
            True,
            ["nF", "nS", "nM", "nG", "nN", "nL", "nE"],
        ),
        (
            "composition",
            False,
            ["nF", "nS", "nM", "nG", "nN"],
        ),
    ],
)
def test_available_meta_properties(mode, sia, expected):
    assert mp.available_meta_properties(mode, sia) == expected


def test_get_mp_object(register_some_mp):
    assert mp.get_meta_property("some_mp").name == "some_mp"


# !!! All tests on meta-properties should use this function for conciseness. !!!
def _test_meta_property(meta_property, string, expected, make_glycan):
    """Test a meta-property on a glycan string.

    Args:
        meta_property (mp.MetaProperty): Meta property function to test.
        string (str): Glycan string. Either GlycoCT or composition.
        expected (Any): Expected result.
        make_glycan (Callable[[str], glyc.Structure]): Factory fixture to make a
            glycan from a string. Either `make_structure` or `make_composition`.
            See `conftest.py` for these two fixtures.
    """
    glycan = make_glycan(string)
    assert meta_property(glycan) == expected


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
    ],
)
def test_count_antenna_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.count_antenna_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 3),
    ],
)
def test_count_fuc_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_fuc_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 1),
        ("H5N4", 0),
        ("H5N4F2E1", 2),
    ],
)
def test_count_fuc_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.count_fuc_mp, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 1),
        (test_glycoct_9, 1),
    ],
)
def test_count_core_fuc_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.count_core_fuc_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 0),
        (test_glycoct_2, 0),
        (test_glycoct_9, 2),
    ],
)
def test_count_antennary_fuc_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.count_antennary_fuc_mp, glycoct, expected, make_structure)


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
    _test_meta_property(mp.count_sia_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H5N4", 0),
        ("H5N4F2E1", 1),
    ],
)
def test_count_sia_comp_mp(comp, expected, make_composition):
    _test_meta_property(mp.count_sia_mp, comp, expected, make_composition)


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
    ],
)
def test_glycan_type_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.glycan_type_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_2, True),
        (test_glycoct_3, False),
        (test_glycoct_4, False),
        (test_glycoct_8, True),
    ],
)
def test_bisection_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.bisection_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 3),
        (test_glycoct_3, 5),
        (test_glycoct_5, 3),
    ],
)
def test_count_man_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_man_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 3),
        ("H6N5", 3),
        ("H5N2", 5),
    ],
)
def test_count_man_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.count_man_mp, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 2),
        (test_glycoct_3, 0),
        (test_glycoct_4, 1),
    ],
)
def test_count_gal_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_gal_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 2),
        ("H6N5", 3),
        ("H5N2", 0),
    ],
)
def test_count_gal_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.count_gal_mp, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, 4),
        (test_glycoct_2, 5),
        (test_glycoct_3, 2),
    ],
)
def test_count_glcnac_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_glcnac_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "comp, expected",
    [
        ("H5N4F1S2", 4),
        ("H6N5", 5),
        ("H5N2", 2),
    ],
)
def test_count_glcnac_mp_comp(comp, expected, make_composition):
    _test_meta_property(mp.count_glcnac_mp, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_1, False),
        (test_glycoct_14, True),
    ],
)
def test_has_poly_lacnac_mp(glycoct, expected, make_structure):
    _test_meta_property(mp.has_poly_lacnac_mp, glycoct, expected, make_structure)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 0),
        (test_glycoct_11, 2),
        (test_glycoct_12, 1),
    ],
)
def test_count_a23_sia_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_a23_sia_mp, glycoct, expected, make_structure)


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
    _test_meta_property(mp.count_a23_sia_mp, comp, expected, make_composition)


@pytest.mark.parametrize(
    "glycoct, expected",
    [
        (test_glycoct_10, 2),
        (test_glycoct_11, 0),
        (test_glycoct_12, 1),
    ],
)
def test_count_a26_sia_mp_struc(glycoct, expected, make_structure):
    _test_meta_property(mp.count_a26_sia_mp, glycoct, expected, make_structure)


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
    _test_meta_property(mp.count_a26_sia_mp, comp, expected, make_composition)
