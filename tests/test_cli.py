from unittest.mock import Mock, MagicMock

import pytest
from click.testing import CliRunner

from glytrait import cli


@pytest.mark.parametrize("sia_linkage", [True, False])
def test_run_workflow(mocker, sia_linkage):
    abund_df_mock = Mock()
    abund_df_mock.columns = ["col1", "col2"]
    formulas_mock = [MagicMock(), MagicMock(), MagicMock()]
    formulas_mock[0].sia_linkage = False
    formulas_mock[1].sia_linkage = True
    formulas_mock[2].sia_linkage = False

    read_input_mock = mocker.patch(
        "glytrait.cli.read_input", return_value=("glycans", abund_df_mock)
    )
    load_formulas_mock = mocker.patch(
        "glytrait.cli.load_formulas", return_value=formulas_mock
    )
    build_meta_property_table_mock = mocker.patch(
        "glytrait.cli.build_meta_property_table", return_value="meta_prop_df"
    )
    calcu_trait_mock = mocker.patch(
        "glytrait.cli.calcu_trait", return_value="trait_df"
    )
    write_output_mock = mocker.patch("glytrait.cli.write_output")

    cli.run_workflow("input_file", "output_file", sia_linkage, "user_formula_file")

    read_input_mock.assert_called_once_with("input_file")
    load_formulas_mock.assert_called_once_with("user_formula_file")
    build_meta_property_table_mock.assert_called_once_with(
        ["col1", "col2"], "glycans", sia_linkage
    )

    if sia_linkage is False:
        formulas_mock.pop(1)
    calcu_trait_mock.assert_called_once_with(
        abund_df_mock, "meta_prop_df", formulas_mock
    )
    write_output_mock.assert_called_once_with(
        "output_file", "trait_df", abund_df_mock, "meta_prop_df", formulas_mock
    )


def test_cli_save_template(clean_dir):
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--save_template", str(clean_dir)])
    assert result.exit_code == 0
    assert "Template saved to" in result.output
    file = clean_dir / "trait_formula.txt"
    content = file.read_text()
    assert "Trait Formula Overview" in content


@pytest.fixture
def empty_file(clean_dir) -> str:
    file = clean_dir / "file"
    file.write_text("")
    return str(file)


def test_cli(mocker, empty_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(empty_file), "output_file"])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(empty_file), "output_file", False, None
    )


def test_cli_sia_linkage(mocker, empty_file):
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(empty_file), "output_file", "-s"])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(empty_file), "output_file", True, None
    )


def test_cli_user_traits(mocker, clean_dir):
    input_file = clean_dir / "file1"
    user_file = clean_dir / "file2"
    input_file.write_text("")
    user_file.write_text("")
    runner = CliRunner()
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow")
    result = runner.invoke(cli.cli, [str(input_file), "output_file", "-f", str(user_file)])
    assert result.exit_code == 0
    run_workflow_mock.assert_called_once_with(
        str(input_file), "output_file", False, str(user_file)
    )


def test_cli_glytrait_error(mocker, empty_file):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow", side_effect=cli.GlyTraitError)
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(empty_file), "output_file"])
    run_workflow_mock.assert_called_once()
    assert result.exit_code != 0
    assert '👎' in result.output


def test_bad_param_error(mocker, empty_file):
    run_workflow_mock = mocker.patch("glytrait.cli.run_workflow", side_effect=cli.GlyTraitError)
    runner = CliRunner()
    result = runner.invoke(cli.cli, [str(empty_file)])
    run_workflow_mock.assert_not_called()
    assert result.exit_code != 0
    assert "You must provide both an input file and an output file." in result.output
