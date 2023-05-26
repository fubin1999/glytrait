import pytest

from glytrait import glycan as glyc


@pytest.fixture
def make_glycan():
    def _make_glycan(string, format="glycoct"):
        return glyc.NGlycan.from_string(string, format=format)
    return _make_glycan


# H5N4S2 (complex)
r"""
Neu5Ac - Gal - Glc2NAc - Man
                            \
                             Man - Glc2NAc - Glc2NAc
                            /
Neu5Ac - Gal - Glc2NAc - Man
"""
test_glycoct_1 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dglc-HEX-1:5
6b:x-dgal-HEX-1:5
7b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
8s:n-acetyl
9s:n-acetyl
10b:x-dman-HEX-1:5
11b:x-dglc-HEX-1:5
12b:x-dgal-HEX-1:5
13b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
14s:n-acetyl
15s:n-acetyl
16s:n-acetyl
17s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:4o(-1+1)5d
5:5o(-1+1)6d
6:6o(-1+2)7d
7:7d(5+1)8n
8:5d(2+1)9n
9:3o(-1+1)10d
10:10o(-1+1)11d
11:11o(-1+1)12d
12:12o(-1+2)13d
13:13d(5+1)14n
14:11d(2+1)15n
15:2d(2+1)16n
16:1d(2+1)17n
"""

# H5N5F1S1 (complex, bisecting)
r"""
Neu5Ac - Gal - Glc2NAc - Man                   Fuc
                            \                   |
                   Glc2NAc - Man - Glc2NAc - Glc2NAc
                            /
         Gal - Glc2NAc - Man
"""
test_glycoct_2 = """RES
1b:x-dglc-HEX-1:5
2b:x-lgal-HEX-1:5|6:d
3b:x-dglc-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dglc-HEX-1:5
6s:n-acetyl
7b:x-dman-HEX-1:5
8b:x-dglc-HEX-1:5
9b:x-dgal-HEX-1:5
10s:n-acetyl
11b:x-dman-HEX-1:5
12b:x-dglc-HEX-1:5
13b:x-dgal-HEX-1:5
14b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
15s:n-acetyl
16s:n-acetyl
17s:n-acetyl
18s:n-acetyl
LIN
1:1o(-1+1)2d
2:1o(-1+1)3d
3:3o(-1+1)4d
4:4o(-1+1)5d
5:5d(2+1)6n
6:4o(-1+1)7d
7:7o(-1+1)8d
8:8o(-1+1)9d
9:8d(2+1)10n
10:4o(-1+1)11d
11:11o(-1+1)12d
12:12o(-1+1)13d
13:13o(-1+2)14d
14:14d(5+1)15n
15:12d(2+1)16n
16:3d(2+1)17n
17:1d(2+1)18n
"""

# H5N2 (high mannose)
r"""
Man
   \
    Man - Man
   /         \
Man           Man - Glc2NAc - Glc2NAc
             /
          Man
"""
test_glycoct_3 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dman-HEX-1:5
6b:x-dman-HEX-1:5
7b:x-dman-HEX-1:5
8s:n-acetyl
9s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:3o(-1+1)5d
5:5o(-1+1)6d
6:5o(-1+1)7d
7:2d(2+1)8n
8:1d(2+1)9n
"""

# H6N3S1 (hybrid)
r"""
               Man
                  \
                   Man - Man
                  /         \
               Man           Man - Glc2NAc - Glc2NAc
                            /
Neu5Ac - Gal - Glc2NAc - Man
"""
test_glycoct_4 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dman-HEX-1:5
6b:x-dman-HEX-1:5
7b:x-dman-HEX-1:5
8b:x-dglc-HEX-1:5
9b:x-dgal-HEX-1:5
10b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
11s:n-acetyl
12s:n-acetyl
13s:n-acetyl
14s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:4o(-1+1)5d
5:4o(-1+1)6d
6:3o(-1+1)7d
7:7o(-1+1)8d
8:8o(-1+1)9d
9:9o(-1+2)10d
10:10d(5+1)11n
11:8d(2+1)12n
12:2d(2+1)13n
13:1d(2+1)14n
"""

# H3N2 (complex)
r"""
Man
   \
    Man - Glc2NAc - Glc2NAc
   /
Man
"""
test_glycoct_5 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dman-HEX-1:5
6s:n-acetyl
7s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:3o(-1+1)5d
5:2d(2+1)6n
6:1d(2+1)7n
"""

# H4N3S1 (complex, monoantennary)
r"""
Neu5Ac - Gal - Glc2NAc - Man
                            \
                             Man - Glc2NAc - Glc2NAc
                            /
                         Man
"""
test_glycoct_6 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dman-HEX-1:5
6b:x-dglc-HEX-1:5
7b:x-dgal-HEX-1:5
8b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
9s:n-acetyl
10s:n-acetyl
11s:n-acetyl
12s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:3o(-1+1)5d
5:5o(-1+1)6d
6:6o(-1+1)7d
7:7o(-1+2)8d
8:8d(5+1)9n
9:6d(2+1)10n
10:2d(2+1)11n
11:1d(2+1)12n
"""

# H6N2 (high mannose)
r"""
Man
   \
    Man - Man
   /         \
Man           Man - Glc2NAc - Glc2NAc
             /
    Man - Man
"""
test_glycoct_7 = """RES
1b:x-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
3b:x-dman-HEX-1:5
4b:x-dman-HEX-1:5
5b:x-dman-HEX-1:5
6b:x-dman-HEX-1:5
7b:x-dman-HEX-1:5
8b:x-dman-HEX-1:5
9s:n-acetyl
10s:n-acetyl
LIN
1:1o(-1+1)2d
2:2o(-1+1)3d
3:3o(-1+1)4d
4:4o(-1+1)5d
5:3o(-1+1)6d
6:6o(-1+1)7d
7:6o(-1+1)8d
8:2d(2+1)9n
9:1d(2+1)10n
"""


class TestGlycan:

    def test_from_string(self):
        glycan = glyc.NGlycan.from_string(test_glycoct_1, format="glycoct")
        assert repr(glycan) == "NGlycan({Gal:2; Man:3; Glc2NAc:4; Neu5Ac:2})"

    def test_from_string_wrong_format(self):
        with pytest.raises(ValueError) as excinfo:
            glyc.NGlycan.from_string(test_glycoct_1, format="unknown")
            assert "Unknown format: unknown" in str(excinfo.value)

    def test_from_string_wrong_string(self):
        with pytest.raises(glyc.StructureParseError) as excinfo:
            glyc.NGlycan.from_string("wrong string", format="glycoct")
            assert "Could not parse string: wrong string" in str(excinfo.value)

    def test_from_glycoct(self):
        glycan = glyc.NGlycan.from_glycoct(test_glycoct_1)
        assert repr(glycan) == "NGlycan({Gal:2; Man:3; Glc2NAc:4; Neu5Ac:2})"

    def test_from_glycoct_wrong_string(self):
        with pytest.raises(glyc.StructureParseError) as excinfo:
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
        ]
    )
    def test_type(self, glycoct, expected, id, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.type == expected

    @pytest.mark.parametrize(
        "glycoct, expected",
        [
            (test_glycoct_1, False),
            (test_glycoct_2, True),
        ]
    )
    def test_bisection(self, glycoct, expected, make_glycan):
        glycan = make_glycan(glycoct)
        assert glycan.is_bisecting == expected

    def test_is_complex(self, make_glycan):
        glycan = make_glycan(test_glycoct_1)
        assert glycan.is_complex() is True

    def test_is_hybrid(self, make_glycan):
        glycan = make_glycan(test_glycoct_4)
        assert glycan.is_hybrid() is True

    def test_is_high_mannose(self, make_glycan):
        glycan = make_glycan(test_glycoct_3)
        assert glycan.is_high_mannose() is True
