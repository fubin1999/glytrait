import pytest

from glytrait import formula as fml
from glytrait import post_filtering as pf


@pytest.mark.parametrize(
    "trait1, trait2, expected",
    [
        ("A2G", "CG", False),
        ("A2Fa", "CFa", True),
        ("A2Fc", "CFc", True),
        ("A2S", "CS", False),
        ("A2SG", "A2G", True),
        ("A2FSG", "A2FG", True),
        ("A2FSG", "A2SG", True),
        ("A2F0G", "A2FG", False),
        ("A2FSG", "A2G", False),
    ],
)
def test_is_child_of(trait1, trait2, expected):
    formulas = fml.load_default_formulas("structure", sia_linkage=True)
    formula_map = {f.name: f for f in formulas}
    formulas1 = formula_map[trait1]
    formulas2 = formula_map[trait2]
    assert pf._is_child_of(formulas1, formulas2) == expected
