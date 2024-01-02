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


@pytest.fixture(autouse=True)
def patch_glytrait(mocker, glytrait_mock):
    mocker.patch("glytrait.cli.GlyTrait", return_value=glytrait_mock)


class TestCLI:
    def test_no_options(self, abundance_file, glycan_file, glytrait_mock):
        runner = CliRunner()
        args = [abundance_file, glycan_file]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called_once_with(
            output_dir=abundance_file.replace(".csv", "_glytrait"),
            abundance_file=abundance_file,
            glycan_file=glycan_file,
            group_file=None,
        )

    @pytest.mark.parametrize("option", ["-g", "--group-file"])
    def test_group_file(
        self, abundance_file, glycan_file, group_file, glytrait_mock, option
    ):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-g", group_file]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once()
        glytrait_mock.run.assert_called_once_with(
            output_dir=abundance_file.replace(".csv", "_glytrait"),
            abundance_file=abundance_file,
            glycan_file=glycan_file,
            group_file=group_file,
        )

    @pytest.mark.parametrize("option", ["-o", "--otuput-dir"])
    def test_output_dir(
        self, abundance_file, glycan_file, output_dir, glytrait_mock, option
    ):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-o", output_dir]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once()
        glytrait_mock.run.assert_called_once_with(
            output_dir=output_dir,
            abundance_file=abundance_file,
            glycan_file=glycan_file,
            group_file=None,
        )

    @pytest.mark.parametrize("option", ["-m", "--mode"])
    def test_mode(self, abundance_file, glycan_file, glytrait_mock, option):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-m", "composition"]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="composition",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-r", "--filter-glycan-ratio"])
    def test_filter_glycan_ratio(
        self, abundance_file, glycan_file, glytrait_mock, option
    ):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-r", "0.5"]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=0.5,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-i", "--impute-method"])
    def test_impute_method(self, abundance_file, glycan_file, glytrait_mock, option):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-i", "mean"]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="mean",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-l", "--sia-linkage"])
    def test_sia_linkage(self, abundance_file, glycan_file, glytrait_mock, option):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-l"]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=True,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-f", "--formula-file"])
    def test_formula_file(
        self, abundance_file, glycan_file, formula_file, glytrait_mock, option
    ):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-f", formula_file]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=formula_file,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize(
        "option, filter",
        [
            ("--filter", True),
            ("--no-filter", False),
        ],
    )
    def test_filter(self, abundance_file, glycan_file, glytrait_mock, option, filter):
        runner = CliRunner()
        args = [abundance_file, glycan_file, option]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=filter,
            correlation_threshold=1.0,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-c", "--corr-threshold"])
    def test_corr_threshold(self, abundance_file, glycan_file, glytrait_mock, option):
        runner = CliRunner()
        args = [abundance_file, glycan_file, "-c", "0.5"]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_called_once_with(
            mode="structure",
            filter_max_na=1.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=0.5,
            sia_linkage=False,
            custom_formula_file=None,
        )
        glytrait_mock.run.assert_called()

    @pytest.mark.parametrize("option", ["-b", "--builtin-formulas"])
    def test_builtin_formulas(self, option, clean_dir, mocker):
        mocker.patch("glytrait.cli.save_builtin_formula", autospec=True)
        dirpath = clean_dir / "formulas"
        dirpath.mkdir()

        runner = CliRunner()
        args = [option, str(dirpath)]
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0
        cli.GlyTrait.assert_not_called()
        cli.save_builtin_formula.assert_called_once_with(str(dirpath))
