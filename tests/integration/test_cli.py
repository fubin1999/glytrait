from pathlib import Path

from click.testing import CliRunner

from glytrait.cli import cli


def copy_data_file(from_: str | Path, to_: str | Path):
    """Copy data file to a new location."""
    with open(from_, "r") as f:
        data = f.read()
    with open(to_, "w") as f:
        f.write(data)


def save_abundance_table(path: str | Path):
    """Save abundance table to a file."""
    data_file = "tests/integration/data/abundance.csv"
    copy_data_file(data_file, path)


def save_structures(path: str | Path):
    """Save structures to a file."""
    data_file = "tests/integration/data/structures.csv"
    copy_data_file(data_file, path)


def save_groups(path: str | Path):
    """Save groups to a file."""
    data_file = "tests/integration/data/groups.csv"
    copy_data_file(data_file, path)


class TestCLI:

    def test_basic(self, tmp_path):
        runner = CliRunner()
        abund_path = str(tmp_path / "abundance.csv")
        struct_path = str(tmp_path / "structures.csv")
        save_abundance_table(abund_path)
        save_structures(struct_path)
        result = runner.invoke(cli, [abund_path, struct_path])
        assert result.exit_code == 0
        assert (tmp_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (tmp_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (tmp_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (tmp_path / "abundance_glytrait/meta_properties.csv").exists()
