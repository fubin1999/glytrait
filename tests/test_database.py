import pytest

from glytrait.database import load_default


@pytest.mark.parametrize("db", ["serum", "IgG"])
def test_load_default(db):
    result = load_default(db, ["H3N3", "H5N4S2"])
    assert len(result) == 2
