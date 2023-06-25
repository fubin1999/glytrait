from pathlib import Path

import pytest
from click.testing import CliRunner

from glytrait import cli
from glytrait import config as gt_config


def test_cli_save_template(clean_dir):
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--save_template", str(clean_dir)])
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


def test_cli(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file)])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == default_config


def test_cli_output_path(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-o", "output_path.xlsx"])
    assert result.exit_code == 0
    config = default_config | dict(output_file="output_path.xlsx")
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_sia_linkage(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "--sia_linkage"])
    assert result.exit_code == 0
    config = default_config | dict(sia_linkage=True)
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_user_traits(mocker, input_file, clean_dir, default_config):
    user_file = clean_dir / "user_formula.txt"
    user_file.write_text("")
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-f", str(user_file)])
    assert result.exit_code == 0
    config = default_config | dict(formula_file=str(user_file))
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_glytrait_error(mocker, input_file):
    run_workflow_mock = mocker.patch(
        "glytrait.cli.run_workflow", side_effect=cli.GlyTraitError
    )
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(input_file)])
    run_workflow_mock.assert_called_once()
    assert result.exit_code != 0
    assert "👎" in result.output


def test_cli_input_file_not_exist(mocker, clean_dir):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    runner = CliRunner()
    input_file = clean_dir / "file.csv"
    result = runner.invoke(cli.cli, [str(input_file)])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_cli_output_dir_not_exist(mocker, input_file, clean_dir, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    output_file = clean_dir / "output_dir" / "output.xlsx"
    result = runner.invoke(cli.cli, [str(input_file), "-o", str(output_file)])
    assert result.exit_code == 0
    config = default_config | dict(output_file=str(output_file))
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_formula_file_not_exist(mocker, input_file, clean_dir):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    formula_file = clean_dir / "formula.txt"
    result = runner.invoke(cli.cli, [str(input_file), "-f", str(formula_file)])
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_cli_filter_off(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "--no-filter"])
    assert result.exit_code == 0
    config = default_config | dict(filter_invalid_traits=False)
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_group_file(mocker, input_file, clean_dir, default_config):
    group_file = clean_dir / "group_file.csv"
    group_file.write_text("")
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-g", str(group_file)])
    assert result.exit_code == 0
    config = default_config | dict(group_file=str(group_file))
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_filter_glycans(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-r", 0.2])
    assert result.exit_code == 0
    config = default_config | dict(filter_glycan_max_na=0.2)
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_impute(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-i", "median"])
    assert result.exit_code == 0
    config = default_config | dict(impute_method="median")
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_structure_file(mocker, input_file, clean_dir, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    structure_file = clean_dir / "structure_file.csv"
    structure_file.touch()
    result = runner.invoke(cli.cli, [str(input_file), "-s", str(structure_file)])
    assert result.exit_code == 0
    config = default_config | dict(structure_file=str(structure_file))
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


def test_cli_databaes(mocker, input_file, default_config):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-d", "serum"])
    assert result.exit_code == 0
    config = default_config | dict(database="serum")
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config


@pytest.mark.parametrize(
    "mode_command, mode",
    [
        ("structure", "structure"),
        ("S", "structure"),
        ("composition", "composition"),
        ("C", "composition"),
    ]
)
def test_cli_mode_structure(mocker, input_file, default_config, mode_command, mode):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-m", mode_command])
    assert result.exit_code == 0
    config = default_config | dict(mode=mode)
    run_workflow_mock.assert_called_once()
    assert run_workflow_mock.call_args[0][0].asdict() == config
