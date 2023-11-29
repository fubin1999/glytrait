import pytest

from glytrait import glycan as glyc
from glytrait.exception import *
from .glycoct import *


class TestGlycan:
    def test_from_string(self):
        glycan = glyc.Structure.from_string("glycan1", test_glycoct_1, format="glycoct")
        assert repr(glycan) == "Structure(name='glycan1')"

    def test_from_string_wrong_format(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.Structure.from_string("glycan1", test_glycoct_1, format="unknown")
            assert "Unknown format: unknown" in str(excinfo.value)

    def test_from_string_wrong_string(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.Structure.from_string("glycan1", "wrong string", format="glycoct")
            assert "Could not parse string: wrong string" in str(excinfo.value)

    def test_from_glycoct(self):
        glycan = glyc.Structure.from_glycoct("glycan1", test_glycoct_1)
        assert repr(glycan) == "Structure(name='glycan1')"

    def test_from_glycoct_wrong_string(self):
        with pytest.raises(StructureParseError) as excinfo:
            glyc.Structure.from_glycoct("glycan1", "wrong string")
            assert "Could not parse string: wrong string" in str(excinfo.value)

    @pytest.mark.parametrize(
        "glycoct, mode, expected",
        [
            (
                test_glycoct_1,
                "bfs",
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
            (
                test_glycoct_1,
                "dfs",
                [
                    "Glc2NAc",
                    "Glc2NAc",
                    "Man",
                    "Man",
                    "Glc2NAc",
                    "Gal",
                    "Neu5Ac",
                    "Man",
                    "Glc2NAc",
                    "Gal",
                    "Neu5Ac",
                ],
            ),
        ],
    )
    def test_traversal(self, glycoct, mode, expected, make_structure):
        glycan = make_structure(glycoct)
        assert [glyc.get_mono_str(mono) for mono in glycan._traversal(mode)] == expected

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
    def test_traversal_skip(self, glycoct, expected, make_structure):
        glycan = make_structure(glycoct)
        it = glycan._traversal("bfs", skip=["Gal", "Neu5Ac"])
        assert list(glyc.get_mono_str(mono) for mono in it) == expected

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
    def test_traversal_only(self, glycoct, expected, make_structure):
        glycan = make_structure(glycoct)
        it = glycan._traversal("bfs", only=["Glc2NAc", "Man"])
        assert list(glyc.get_mono_str(mono) for mono in it) == expected

    def test_traversal_both_skip_and_only(self, make_structure):
        glycan = make_structure(test_glycoct_1)
        with pytest.raises(ValueError) as excinfo:
            for _ in glycan._traversal("bfs", skip=["Gal"], only=["Glc2NAc"]):
                pass
        assert "Cannot specify both `skip` and `only`." in str(excinfo.value)

    def test_wrong_traversal_mode(self, make_structure):
        glycan = make_structure(test_glycoct_1)
        with pytest.raises(ValueError) as excinfo:
            for _ in glycan._traversal("wrong"):
                pass
        assert "Unknown traversal method: wrong" in str(excinfo.value)

    def test_get(self, make_structure):
        glycan = make_structure(test_glycoct_1)
        assert glycan.get("Glc2NAc") == 4
        assert glycan.get("Man") == 3


class TestComposition:
    def test_init_unknown_mono(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition("H5N4P1", dict(H=5, N=4, P=1))
        assert "Unknown monosaccharide: P" in str(excinfo.value)

    def test_init_negative_mono(self):
        with pytest.raises(CompositionParseError) as excinfo:
            glyc.Composition("H5N4F-1", dict(H=5, N=4, F=-1))
        assert "Monosacharride must be above 0: F=-1." in str(excinfo.value)

    @pytest.mark.parametrize(
        "comp, expected",
        [
            ("H5N4F1S1", {"H": 5, "N": 4, "F": 1, "S": 1}),
            ("H5N4F1", {"H": 5, "N": 4, "F": 1}),
            ("N4F1S1H5", {"H": 5, "N": 4, "F": 1, "S": 1}),
            ("H5N4F1E1L1", {"H": 5, "N": 4, "F": 1, "E": 1, "L": 1}),
            ("H5N4F1E1", {"H": 5, "N": 4, "F": 1, "E": 1}),
        ],
    )
    def test_from_string(self, comp, expected):
        result = glyc.Composition.from_string(comp)
        assert result.asdict() == expected

    def test_from_string_name(self):
        result = glyc.Composition.from_string("H5N4F1S1")
        assert result.name == "H5N4F1S1"

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

    def test_get(self):
        comp = glyc.Composition("H5N4F1S1", dict(H=5, N=4, F=1, S=1))
        assert comp.get("H") == 5
        assert comp.get("F") == 1
        assert comp.get("S") == 1
        assert comp.get("N") == 4


def test_load_glycans():
    names = ["test_1", "test_2", "test_3"]
    glycocts = [test_glycoct_1, test_glycoct_2, test_glycoct_3]
    glycans = glyc.load_structures(names, glycocts)
    assert len(glycans) == 3


def test_load_compositions_no_sia_linkage():
    comps = ["H5N4", "H5N2", "H5N4S1"]
    result = glyc.load_compositions(comps)
    assert len(result) == 3


def test_load_compositions_sia_linkage():
    comps = ["H5N4L1", "H5N4E1", "H5N4E1L1"]
    result = glyc.load_compositions(comps)
    assert len(result) == 3
