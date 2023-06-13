from pathlib import Path

import pytest
from click.testing import CliRunner

from glytrait import cli


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


def test_cli(mocker, input_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file)])
    output_path = input_file.with_name(input_file.stem + "_glytrait.xlsx")
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_path), False, None, True, None
    )


def test_cli_output_path(mocker, input_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-o", "output_path.xlsx"])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), "output_path.xlsx", False, None, True, None
    )


def test_cli_sia_linkage(mocker, input_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-s"])
    output_path = input_file.with_name(input_file.stem + "_glytrait.xlsx")
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_path), True, None, True, None
    )


def test_cli_user_traits(mocker, input_file, clean_dir):
    user_file = clean_dir / "user_formula.txt"
    user_file.write_text("")
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-f", str(user_file)])
    output_path = input_file.with_name(input_file.stem + "_glytrait.xlsx")
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_path), False, str(user_file), True, None
    )


def test_cli_glytrait_error(mocker, input_file):
    run_workflow_mock = mocker.patch(
        "glytrait.cli.run_workflow", side_effect=cli.GlyTraitError
    )
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(input_file)])
    run_workflow_mock.assert_called_once()
    assert result.exit_code != 0
    assert "👎" in result.output


def test_cli_input_not_csv(mocker, clean_dir):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    runner = CliRunner()
    input_file = clean_dir / "file.not_csv"
    input_file.write_text("")
    result = runner.invoke(cli.cli, [str(input_file)])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "Input file must be a .csv file." in result.output


def test_cli_output_not_xlsx(mocker, input_file):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(input_file), "-o", "output_path.not_xlsx"])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "Output file must be a .xlsx file." in result.output


def test_cli_formula_not_txt(mocker, input_file, clean_dir):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    runner = CliRunner()
    formula_file = clean_dir / "formula.not_txt"
    formula_file.write_text("")
    result = runner.invoke(cli.cli, [str(input_file), "-f", str(formula_file)])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "Formula file must be a .txt file." in result.output


def test_cli_input_file_not_exist(mocker, clean_dir):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    runner = CliRunner()
    input_file = clean_dir / "file.csv"
    result = runner.invoke(cli.cli, [str(input_file)])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_cli_output_dir_not_exist(mocker, input_file, clean_dir):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    output_file = clean_dir / "output_dir" / "output.xlsx"
    result = runner.invoke(cli.cli, [str(input_file), "-o", str(output_file)])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_file), False, None, True, None
    )


def test_cli_formula_file_not_exist(mocker, input_file, clean_dir):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    formula_file = clean_dir / "formula.txt"
    result = runner.invoke(cli.cli, [str(input_file), "-f", str(formula_file)])
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_cli_filter_off(mocker, input_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "--no-filter"])
    output_path = input_file.with_name(input_file.stem + "_glytrait.xlsx")
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_path), False, None, False, None
    )


def test_cli_group_file(mocker, input_file, clean_dir):
    group_file = clean_dir / "group_file.csv"
    group_file.write_text("")
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "-g", str(group_file)])
    output_path = input_file.with_name(input_file.stem + "_glytrait.xlsx")
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), str(output_path), False, None, True, str(group_file)
    )
