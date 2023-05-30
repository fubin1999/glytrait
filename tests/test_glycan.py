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
        ]
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
        ]
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
        ]
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
        ]
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
        ]
    )
    def test_count_gal(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.count_gal() == expected

    def test_init_not_Nglycan(self, make_glycan):
        with pytest.raises(StructureParseError) as excinfo:
            make_glycan(test_glycoct_13)
        assert "This is not an N-glycan." in str(excinfo.value)
