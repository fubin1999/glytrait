import pytest


@pytest.fixture
def clean_dir(tmp_path):
    d = tmp_path / "clean_dir"
    d.mkdir()
    return d
