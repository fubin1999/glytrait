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
    file_dir = Path("tests/integration/data")
    copy_data_file(file_dir / "abundance.csv", tmp_path / "abundance.csv")
    copy_data_file(file_dir / "structures.csv", tmp_path / "structures.csv")
    copy_data_file(file_dir / "groups.csv", tmp_path / "groups.csv")
    copy_data_file(file_dir / "formulas.txt", tmp_path / "formulas.txt")
    copy_data_file(file_dir / "compositions.csv", tmp_path / "compositions.csv")
    copy_data_file(file_dir / "meta_properties.csv", tmp_path / "meta_properties.csv")
    return tmp_path


def assert_basic_result_files_exist(result_dir: Path) -> None:
    """Assert 'derived_traits.csv', 'derived_traits_filtered.csv',
    'glycan_abundance_processed.csv', and 'meta_properties.csv' exist."""
    assert (result_dir / "derived_traits.csv").exists()
    assert (result_dir / "derived_traits_filtered.csv").exists()
    assert (result_dir / "glycan_abundance_processed.csv").exists()
    assert (result_dir / "meta_properties.csv").exists()


class TestCLI:

    def test_basic(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)
        df = pd.read_csv(result_dir / "derived_traits.csv", index_col=0)
        assert "CE" not in df.columns

    def test_run_with_mp_file(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--mp-file",
            str(prepared_path / "meta_properties.csv")
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)

    def test_group(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "-g",
            str(prepared_path / "groups.csv"),
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert (result_dir / "anova.csv").exists()
        assert (result_dir / "post_hoc.csv").exists()
        assert_basic_result_files_exist(result_dir)

    def test_formula_file(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "-f",
            str(prepared_path / "formulas.txt"),
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)
        df = pd.read_csv(result_dir / "derived_traits.csv", index_col=0)
        assert df.shape[1] == 8

    def test_mode(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "compositions.csv"),
            "-m",
            "C",
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)
        df = pd.read_csv(result_dir / "derived_traits.csv", index_col=0)
        assert "Hb" in df.columns

    def test_output(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "-o",
            str(prepared_path / "output"),
        ]
        result_dir = prepared_path / "output"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)

    def test_sia_linkage(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "--sia-linkage",
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert_basic_result_files_exist(result_dir)
        df = pd.read_csv(result_dir / "derived_traits.csv", index_col=0)
        assert "CE" in df.columns

    def test_no_filter(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "--no-filter",
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert (result_dir / "derived_traits.csv").exists()
        assert (result_dir / "glycan_abundance_processed.csv").exists()
        assert (result_dir / "meta_properties.csv").exists()

    def test_no_filter_has_groups(self, prepared_path):
        runner = CliRunner()
        args = [
            str(prepared_path / "abundance.csv"),
            "--glycan-file",
            str(prepared_path / "structures.csv"),
            "-g",
            str(prepared_path / "groups.csv"),
            "--no-filter",
        ]
        result_dir = prepared_path / "abundance_glytrait"

        result = runner.invoke(cli, args)

        assert result.exit_code == 0
        assert (result_dir / "derived_traits.csv").exists()
        assert (result_dir / "glycan_abundance_processed.csv").exists()
        assert (result_dir / "meta_properties.csv").exists()
        assert not (result_dir / "anova.csv").exists()
        assert "Warning" in result.output

    def test_built_in_formulas(self, tmp_path):
        runner = CliRunner()
        args = ["-b", str(tmp_path)]
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert (tmp_path / "struc_builtin_formulas.txt").exists()
        assert (tmp_path / "comp_builtin_formulas.txt").exists()

    def test_welcome_msg(self):
        runner = CliRunner()
        result = runner.invoke(cli, [])
        assert result.exit_code == 0
        assert "Welcome to GlyTrait!" in result.output
