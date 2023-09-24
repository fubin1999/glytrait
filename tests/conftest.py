import pytest

from glytrait import glycan as glyc


@pytest.fixture
def clean_dir(tmp_path):
    d = tmp_path / "clean_dir"
    d.mkdir()
    return d


@pytest.fixture
def make_structure():
    def _make_structure(string, format="glycoct"):
        return glyc.Structure.from_string("glycan", string, format=format)
    return _make_structure


@pytest.fixture
def make_composition():
    def _make_composition(string):
        return glyc.Composition.from_string(string)
    return _make_composition
