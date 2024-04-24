import pytest
from click.testing import CliRunner

from glytrait import cli


@pytest.fixture
def abundance_file(clean_dir) -> str:
    file = clean_dir / "abundance.csv"
    file.touch()
    return str(file)


@pytest.fixture
def glycan_file(clean_dir) -> str:
    file = clean_dir / "glycan.csv"
    file.touch()
    return str(file)


@pytest.fixture
def group_file(clean_dir) -> str:
    file = clean_dir / "group.csv"
    file.touch()
    return str(file)


@pytest.fixture
def formula_file(clean_dir) -> str:
    file = clean_dir / "formulas.txt"
    file.touch()
    return str(file)


@pytest.fixture
def output_dir(clean_dir) -> str:
    path = clean_dir / "output"
    path.mkdir()
    return str(path)


@pytest.fixture
def glytrait_mock(mocker):
    glytrait = mocker.Mock()
    glytrait.run = mocker.Mock()
    return glytrait
