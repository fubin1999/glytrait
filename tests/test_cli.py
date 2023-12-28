from pathlib import Path

import pytest
from click.testing import CliRunner

from glytrait import cli, config as gt_config

pytestmark = pytest.mark.skip(reason="CLI is not ready yet.")


def test_cli_save_template(clean_dir):
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--save-template", str(clean_dir)])
    assert result.exit_code == 0
    assert "Template saved to" in result.output
    file = clean_dir / "trait_formula.txt"
    content = file.read_text()
    assert "Trait Formula Overview" in content


@pytest.fixture
def input_file(clean_dir) -> Path:
    """Create a csv file in a clean temperary directory."""
    file = clean_dir / "file.csv"
    file.write_text("")
    return file


@pytest.fixture
def output_file(input_file) -> Path:
    """Create a csv file in a clean temperary directory."""
    return input_file.with_name(input_file.stem + "_glytrait.xlsx")


@pytest.fixture
def default_config(input_file, output_file) -> dict:
    cfg = gt_config.default_config.copy()
    cfg.update(
        input_file=str(input_file),
        output_file=str(output_file),
        mode="structure",
    )
    return cfg


@pytest.fixture
def workflow_class_mock(mocker):
    return mocker.patch("glytrait.cli.Workflow", autospec=True)


@pytest.fixture
def workflow_mock(workflow_class_mock):
    return workflow_class_mock.return_value


def _test_normal_cases(workflow_mock, workflow_class_mock, config, args):
    """Test normal cases when CLI is invaked with proper arguments."""
    runner = CliRunner()
    result = runner.invoke(cli.cli, args)
    assert result.exit_code == 0
    workflow_mock.run.assert_called_once()
    assert workflow_class_mock.call_args[0][0].asdict() == config


def _test_invalid_args(
    workflow_mock,
    workflow_class_mock,
    args,
    msg,
):
    """Test cases when CLI raised with invalid args."""
    runner = CliRunner()
    result = runner.invoke(cli.cli, args)
    assert result.exit_code != 0
    assert msg in result.output
    workflow_mock.run.assert_not_called()
    workflow_class_mock.assert_not_called()


def test_cli(input_file, default_config, workflow_class_mock, workflow_mock):
    args = [str(input_file)]
    _test_normal_cases(workflow_mock, workflow_class_mock, default_config, args)


def test_cli_output_path(
    input_file, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "-o", "output_path.xlsx"]
    config = default_config | dict(output_file="output_path.xlsx")
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_sia_linkage(
    input_file, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "-l"]
    config = default_config | dict(sia_linkage=True)
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_user_traits(
    input_file, clean_dir, default_config, workflow_class_mock, workflow_mock
):
    user_file = clean_dir / "user_formula.txt"
    user_file.write_text("")
    args = [str(input_file), "-f", str(user_file)]
    config = default_config | dict(formula_file=str(user_file))
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_glytrait_error(
    input_file, default_config, workflow_mock, workflow_class_mock
):
    workflow_mock.run.side_effect = cli.GlyTraitError
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(input_file)])
    workflow_mock.run.assert_called_once()
    assert workflow_class_mock.call_args[0][0].asdict() == default_config
    assert result.exit_code != 0
    assert "👎" in result.output


def test_cli_input_file_not_exist(
    clean_dir, default_config, workflow_class_mock, workflow_mock
):
    input_file = clean_dir / "not_exist.csv"
    args = [str(input_file)]
    _test_invalid_args(workflow_mock, workflow_class_mock, args, "does not exist")


def test_cli_output_dir_not_exist(
    input_file, clean_dir, default_config, workflow_class_mock, workflow_mock
):
    output_file = clean_dir / "output_dir" / "output.xlsx"
    args = [str(input_file), "-o", str(output_file)]
    config = default_config | dict(output_file=str(output_file))
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_formula_file_not_exist(
    input_file, clean_dir, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "-f", "not_exist.txt"]
    _test_invalid_args(workflow_mock, workflow_class_mock, args, "does not exist")


def test_cli_filter_off(input_file, default_config, workflow_class_mock, workflow_mock):
    args = [str(input_file), "--no-filter"]
    config = default_config | dict(post_filtering=False)
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_filter_glycans(
    input_file, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "-r", 0.2]
    config = default_config | dict(filter_glycan_max_na=0.2)
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_impute(input_file, default_config, workflow_class_mock, workflow_mock):
    args = [str(input_file), "-i", "median"]
    config = default_config | dict(impute_method="median")
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_structure_file(
    input_file, clean_dir, default_config, workflow_class_mock, workflow_mock
):
    structure_file = clean_dir / "structure_file.csv"
    structure_file.touch()
    args = [str(input_file), "-s", str(structure_file)]
    config = default_config | dict(structure_file=str(structure_file))
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_databaes(input_file, default_config, workflow_class_mock, workflow_mock):
    args = [str(input_file), "-d", "serum"]
    config = default_config | dict(database="serum")
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


@pytest.mark.parametrize(
    "mode_command, mode",
    [
        ("structure", "structure"),
        ("S", "structure"),
        ("composition", "composition"),
        ("C", "composition"),
    ],
)
def test_cli_mode_structure(
    input_file, default_config, mode_command, mode, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "-m", mode_command]
    config = default_config | dict(mode=mode)
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_save_built_in_formulas(clean_dir):
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["-b", str(clean_dir)])
    assert result.exit_code == 0
    assert "Built-in formulas saved to" in result.output
    file1 = clean_dir / "struc_builtin_formulas.txt"
    file2 = clean_dir / "comp_builtin_formulas.txt"
    assert file1.exists()
    assert file2.exists()


def test_cli_corr_threshold(
    clean_dir, input_file, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "--corr-threshold", 0.5]
    config = default_config | dict(corr_threshold=0.5)
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)


def test_cli_corr_method(
    clean_dir, input_file, default_config, workflow_class_mock, workflow_mock
):
    args = [str(input_file), "--corr-method", "spearman"]
    config = default_config | dict(corr_method="spearman")
    _test_normal_cases(workflow_mock, workflow_class_mock, config, args)
