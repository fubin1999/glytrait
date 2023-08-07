import pytest

from glytrait import glycan as glyc


@pytest.fixture
def clean_dir(tmp_path):
    d = tmp_path / "clean_dir"
    d.mkdir()
    return d


@pytest.fixture
def make_glycan():
    def _make_glycan(string, format="glycoct"):
        return glyc.NGlycan.from_string("glycan", string, format=format)

    return _make_glycan
