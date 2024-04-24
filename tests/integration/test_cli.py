from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from glytrait.cli import cli


def copy_data_file(from_: str | Path, to_: str | Path):
    """Copy data file to a new location."""
    with open(from_, "r") as f:
        data = f.read()
    with open(to_, "w") as f:
        f.write(data)


@pytest.fixture
def prepared_path(tmp_path) -> Path:
    """Prepare data files."""
    copy_data_file("tests/integration/data/abundance.csv", tmp_path / "abundance.csv")
    copy_data_file("tests/integration/data/structures.csv", tmp_path / "structures.csv")
    copy_data_file("tests/integration/data/groups.csv", tmp_path / "groups.csv")
    copy_data_file("tests/integration/data/formulas.txt", tmp_path / "formulas.txt")
    copy_data_file("tests/integration/data/compositions.csv", tmp_path / "compositions.csv")
    return tmp_path


class TestCLI:

    def test_basic(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()
        df = pd.read_csv(prepared_path / "abundance_glytrait/derived_traits.csv", index_col=0)
        assert "CE" not in df.columns

    def test_group(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "-g",
            str(prepared_path / "groups.csv"),
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/anova.csv").exists()
        assert (prepared_path / "abundance_glytrait/post_hoc.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()

    def test_formula_file(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "-f",
            str(prepared_path / "formulas.txt"),
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()
        df = pd.read_csv(prepared_path / "abundance_glytrait/derived_traits.csv", index_col=0)
        assert df.shape[1] == 8

    def test_mode(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "compositions.csv"),
            "-m",
            "C",
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()
        df = pd.read_csv(prepared_path / "abundance_glytrait/derived_traits.csv", index_col=0)
        assert "Hb" in df.columns

    def test_output(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "-o",
            str(prepared_path / "output"),
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "output/derived_traits.csv").exists()
        assert (prepared_path / "output/derived_traits_filtered.csv").exists()
        assert (prepared_path / "output/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "output/meta_properties.csv").exists()

    def test_sia_linkage(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "--sia-linkage",
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/derived_traits_filtered.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()
        df = pd.read_csv(prepared_path / "abundance_glytrait/derived_traits.csv", index_col=0)
        assert "CE" in df.columns

    def test_no_filter(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "--no-filter",
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()

    def test_no_filter_has_groups(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            str(prepared_path / "structures.csv"),
            "-g",
            str(prepared_path / "groups.csv"),
            "--no-filter",
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (prepared_path / "abundance_glytrait/derived_traits.csv").exists()
        assert (prepared_path / "abundance_glytrait/glycan_abundance_processed.csv").exists()
        assert (prepared_path / "abundance_glytrait/meta_properties.csv").exists()
        assert not (prepared_path / "abundance_glytrait/anova.csv").exists()
        assert "Warning" in result.output

    def test_built_in_formulas(self, tmp_path):
        runner = CliRunner()
        args = [
            "-b",
            str(tmp_path)
        ]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (tmp_path / "struc_builtin_formulas.txt").exists()
        assert (tmp_path / "comp_builtin_formulas.txt").exists()

    def test_welcome_msg(self):
        runner = CliRunner()
        result = runner.invoke(cli, [])
        assert result.exit_code == 0
        assert "Welcome to GlyTrait!" in result.output
