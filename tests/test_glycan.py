import pytest

from glytrait import glycan as glyc
from glytrait.exception import *
from .glycoct import *


class TestGlycan:
    def test_from_string(self):
        glycan = glyc.NGlycan.from_string(test_glycoct_1, format="glycoct")
        assert repr(glycan) == "NGlycan({Gal:2; Man:3; Glc2NAc:4; Neu5Ac:2})"

    def test_from_string_wrong_format(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.NGlycan.from_string(test_glycoct_1, format="unknown")
            assert "Unknown format: unknown" in str(excinfo.value)

    def test_from_string_wrong_string(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.NGlycan.from_string("wrong string", format="glycoct")
            assert "Could not parse string: wrong string" in str(excinfo.value)

    def test_from_glycoct(self):
        glycan = glyc.NGlycan.from_glycoct(test_glycoct_1)
        assert repr(glycan) == "NGlycan({Gal:2; Man:3; Glc2NAc:4; Neu5Ac:2})"

    def test_from_glycoct_wrong_string(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.NGlycan.from_glycoct("wrong string")
            assert "Could not parse string: wrong string" in str(excinfo.value)

    @pytest.mark.parametrize(
        "id, glycoct, expected",
        [
            (1, test_glycoct_1, glyc.GlycanType.COMPLEX),
            (2, test_glycoct_2, glyc.GlycanType.COMPLEX),
            (3, test_glycoct_3, glyc.GlycanType.HIGH_MANNOSE),
            (4, test_glycoct_4, glyc.GlycanType.HYBRID),
            (5, test_glycoct_5, glyc.GlycanType.COMPLEX),
            (6, test_glycoct_6, glyc.GlycanType.COMPLEX),
            (7, test_glycoct_7, glyc.GlycanType.HIGH_MANNOSE),
        ],
    )
    def test_type(self, glycoct, expected, id, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.type == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (
                test_glycoct_1,
                [
                    "Glc2NAc",
                    "Glc2NAc",
                    "Man",
                    "Man",
                    "Man",
                    "Glc2NAc",
                    "Glc2NAc",
                    "Gal",
                    "Gal",
                    "Neu5Ac",
                    "Neu5Ac",
                ],
            ),
        ],
    )
    def test_traversal(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert (
            list(glyc.get_mono_comp(mono) for mono in glycan._traversal("bfs"))
            == expected
        )

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (
                test_glycoct_1,
                [
                    "Glc2NAc",
                    "Glc2NAc",
                    "Man",
                    "Man",
                    "Man",
                    "Glc2NAc",
                    "Glc2NAc",
                ],
            )
        ],
    )
    def test_traversal_skip(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        it = glycan._traversal("bfs", skip=["Gal", "Neu5Ac"])
        assert list(glyc.get_mono_comp(mono) for mono in it) == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (
                test_glycoct_1,
                [
                    "Glc2NAc",
                    "Glc2NAc",
                    "Man",
                    "Man",
                    "Man",
                    "Glc2NAc",
                    "Glc2NAc",
                ],
            )
        ],
    )
    def test_traversal_only(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        it = glycan._traversal("bfs", only=["Glc2NAc", "Man"])
        assert list(glyc.get_mono_comp(mono) for mono in it) == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, False),
            (test_glycoct_2, True),
        ],
    )
    def test_bisection(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.is_bisecting() == expected

    def test_is_complex(self, make_glycan):
        glycan = make_glycan(test_glycoct_1)
        assert glycan.is_complex() is True

    def test_is_hybrid(self, make_glycan):
        glycan = make_glycan(test_glycoct_4)
        assert glycan.is_hybrid() is True

    def test_is_high_mannose(self, make_glycan):
        glycan = make_glycan(test_glycoct_3)
        assert glycan.is_high_mannose() is True

    @pytest.mark.parametrize(
        """glycoct, expected""",
        [
            (test_glycoct_1, 2),
            (test_glycoct_2, 2),
            (test_glycoct_3, 0),
            (test_glycoct_5, 0),
            (test_glycoct_6, 1),
            (test_glycoct_8, 3),
        ],
    )
    def test_count_branches(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_antenna() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 0),
            (test_glycoct_2, 1),
            (test_glycoct_3, 0),
            (test_glycoct_4, 0),
            (test_glycoct_5, 0),
            (test_glycoct_6, 0),
            (test_glycoct_7, 0),
            (test_glycoct_8, 1),
            (test_glycoct_9, 3),
        ],
    )
    def test_count_fuc(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_fuc() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 0),
            (test_glycoct_2, 1),
            (test_glycoct_8, 1),
            (test_glycoct_9, 1),
        ],
    )
    def test_count_core_fuc(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_core_fuc() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 0),
            (test_glycoct_2, 0),
            (test_glycoct_8, 0),
            (test_glycoct_9, 2),
        ],
    )
    def test_count_antennary_fuc(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_antennary_fuc() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 2),
            (test_glycoct_2, 1),
            (test_glycoct_3, 0),
            (test_glycoct_4, 1),
            (test_glycoct_5, 0),
            (test_glycoct_6, 1),
            (test_glycoct_7, 0),
            (test_glycoct_8, 1),
            (test_glycoct_9, 1),
        ],
    )
    def test_count_sia(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_sia() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_10, 0),
            (test_glycoct_11, 2),
            (test_glycoct_12, 1),
            (test_glycoct_3, 0),
        ],
    )
    def test_count_a23_sia(self, glycoct, expected, make_glycan):
        assert make_glycan(glycoct).count_a23_sia() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_10, 2),
            (test_glycoct_11, 0),
            (test_glycoct_12, 1),
            (test_glycoct_3, 0),
        ],
    )
    def test_count_a26_sia(self, glycoct, expected, make_glycan):
        assert make_glycan(glycoct).count_a26_sia() == expected

    def test_count_a23_a26_sia_no_linkage_info(self, make_glycan):
        glycan = make_glycan(test_glycoct_1)
        with pytest.raises(SiaLinkageError):
            glycan.count_a23_sia()
        with pytest.raises(SiaLinkageError):
            glycan.count_a26_sia()

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 3),
            (test_glycoct_3, 5),
            (test_glycoct_4, 5),
        ],
    )
    def test_count_man(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_man() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, 2),
            (test_glycoct_2, 2),
            (test_glycoct_3, 0),
            (test_glycoct_4, 1),
        ],
    )
    def test_count_gal(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_gal() == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, False),
            (test_glycoct_2, False),
            (test_glycoct_8, False),
            (test_glycoct_14, True),
        ],
    )
    def test_has_poly_lacnac(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.has_poly_lacnac() == expected

    def test_init_not_Nglycan(self, make_glycan):
        with pytest.raises(StructureParseError) as excinfo:
            make_glycan(test_glycoct_13)
        assert "This is not an N-glycan." in str(excinfo.value)


class TestComposition:
    def test_init_unknown_mono(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition(dict(H=5, N=4, P=1), sia_linkage=False)
        assert "Unknown monosaccharide: P" in str(excinfo.value)

    def test_init_negative_mono(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition(dict(H=5, N=4, F=-1), sia_linkage=False)
        assert "Monosacharride must be above 0: F=-1." in str(excinfo.value)

    @pytest.mark.parametrize(
        "comp, sia_linkage, expected",
        [
            ("H5N4F1S1", False, {"H": 5, "N": 4, "F": 1, "S": 1}),
            ("H5N4F1", False, {"H": 5, "N": 4, "F": 1, "S": 0}),
            ("N4F1S1H5", False, {"H": 5, "N": 4, "F": 1, "S": 1}),
            ("H5N4F1E1L1", True, {"H": 5, "N": 4, "F": 1, "E": 1, "L": 1}),
            ("H5N4F1E1", True, {"H": 5, "N": 4, "F": 1, "E": 1, "L": 0}),
        ],
    )
    def test_from_string(self, comp, sia_linkage, expected):
        result = glyc.Composition.from_string(comp, sia_linkage=sia_linkage)
        assert result.asdict() == expected

    @pytest.mark.parametrize("s", ["H5N4FS", "1HN4", "abc"])
    def test_from_string_invalid(self, s):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition.from_string(s)
        assert f"Invalid composition: {s}." in str(excinfo.value)

    def test_from_string_invalid_mono(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition.from_string("H5N4F1P1")
        assert "Unknown monosaccharide: P" in str(excinfo.value)

    def test_from_string_emtpy(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition.from_string("")
        assert "Empty string." in str(excinfo.value)

    def test_from_string_invalid_sia_linkage(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition.from_string("H5N4F1S1", sia_linkage=True)
        msg = "'S' is not allow for sialic-acid-linkage-specified composition."
        assert msg in str(excinfo.value)

    def test_from_string_invalid_sia_linkage_2(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition.from_string("H5N4F1E1", sia_linkage=False)
        msg = (
            "'E' and 'L' is not allow for sialic-acid-linkage-unspecified composition."
        )
        assert msg in str(excinfo.value)

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4", False),
            ("H5N3", False),
            ("H5N5", True),
        ],
    )
    def test_is_high_branching(self, comp, expected):
        assert glyc.Composition.from_string(comp).is_high_branching() == expected

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4", True),
            ("H5N3", True),
            ("H5N5", False),
        ],
    )
    def test_is_low_branching(self, comp, expected):
        assert glyc.Composition.from_string(comp).is_low_branching() == expected

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N2", 0),
            ("H4N2", 0),
            ("H4N3", 1),
            ("H5N5", 2),
            ("H6N6", 3),
            ("H7N6", 4),
        ],
    )
    def test_count_gal(self, comp, expected):
        assert glyc.Composition.from_string(comp).count_gal() == expected

    @pytest.mark.parametrize(
        "comp, sia_linkage, expected",
        [
            ("H5N4", False, 0),
            ("H5N4F1", False, 0),
            ("H5N4S1", False, 1),
            ("H5N4S2", False, 2),
            ("H5N4E2", True, 2),
            ("H5N4E1L1", True, 2),
        ],
    )
    def test_count_sia(self, comp, sia_linkage, expected):
        result = glyc.Composition.from_string(comp, sia_linkage=sia_linkage).count_sia()
        assert result == expected

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4", 0),
            ("H5N4F1", 1),
            ("H5N4F2", 2),
            ("H5N4F2S1", 2),
        ],
    )
    def test_count_fuc(self, comp, expected):
        assert glyc.Composition.from_string(comp).count_fuc() == expected

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4L1", 1),
            ("H5N4E1", 0),
            ("H5N4E1L1", 1),
        ],
    )
    def test_count_a23_sia(self, comp, expected):
        result = glyc.Composition.from_string(comp, sia_linkage=True).count_a23_sia()
        assert result == expected

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4L1", 0),
            ("H5N4E1", 1),
            ("H5N4E1L1", 1),
        ],
    )
    def test_count_a26_sia(self, comp, expected):
        result = glyc.Composition.from_string(comp, sia_linkage=True).count_a26_sia()
        assert result == expected


def test_load_glycans():
    glycocts = [test_glycoct_1, test_glycoct_2, test_glycoct_3]
    glycans = glyc.load_glycans(glycocts)
    assert len(glycans) == 3


def test_load_compositions_no_sia_linkage():
    comps = ["H5N4", "H5N2", "H5N4S1"]
    result = glyc.load_compositions(comps, sia_linkage=False)
    assert len(result) == 3


def test_load_compositions_sia_linkage():
    comps = ["H5N4L1", "H5N4E1", "H5N4E1L1"]
    result = glyc.load_compositions(comps, sia_linkage=True)
    assert len(result) == 3
