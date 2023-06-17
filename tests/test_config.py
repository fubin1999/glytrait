import pytest

from glytrait import config


class TestConfig:
    @pytest.fixture
    def input_file(self, clean_dir) -> str:
        input_file = clean_dir / "test.csv"
        input_file.touch()
        return str(input_file)

    @pytest.fixture
    def output_file(self, clean_dir) -> str:
        output_file = clean_dir / "test.xlsx"
        return str(output_file)

    def test_init(self, input_file, output_file):
        new = {"input_file": input_file, "output_file": output_file}
        cfg = config.Config(new)
        assert cfg.asdict() == config.default_config | new

    def test_init_without_input_file(self, output_file):
        new = {"output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Missing config key: input_file" in str(excinfo.value)

    def test_input_file_not_csv(self, clean_dir, output_file):
        input_file = clean_dir / "test.txt"
        input_file.touch()
        input_file = str(input_file)
        new = {"input_file": input_file, "output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a CSV file" in str(excinfo.value)

    def test_input_file_not_exists(self, clean_dir, output_file):
        input_file = str(clean_dir / "test.csv")
        new = {"input_file": input_file, "output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file does not exist" in str(excinfo.value)

    def test_input_not_a_file(self, clean_dir, output_file):
        input_file = clean_dir / "test.csv"
        input_file.mkdir()
        input_file = str(input_file)
        new = {"input_file": input_file, "output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a file, not directory" in str(excinfo.value)

    def test_input_file_not_str(self, output_file):
        new = {"input_file": 1, "output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a string" in str(excinfo.value)

    def test_output_file_not_xlsx(self, clean_dir, input_file):
        output_file = str(clean_dir / "test.txt")
        new = {"input_file": input_file, "output_file": output_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Output file must be a XLSX file" in str(excinfo.value)

    def test_filter_glycan_max_na(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "filter_glycan_max_na": 0.5,
        }
        cfg = config.Config(new)
        assert cfg.get("filter_glycan_max_na") == 0.5

    def test_filter_glycan_max_na_out_of_range(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "filter_glycan_max_na": 1.5,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "filter_glycan_max_na must be between 0 and 1" in str(excinfo.value)

    def test_filter_glycan_max_na_str(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "filter_glycan_max_na": "0.5",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "filter_glycan_max_na must be a float" in str(excinfo.value)

    def test_filter_glycan_max_na_int(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "filter_glycan_max_na": 1,
        }
        cfg = config.Config(new)
        assert cfg.get("filter_glycan_max_na") == 1

    def test_impute_method(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "impute_method": "min",
        }
        cfg = config.Config(new)
        assert cfg.get("impute_method") == "min"

    def test_impute_method_not_supported(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "impute_method": "max",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "impute_method must be one of" in str(excinfo.value)

    def test_impute_method_not_str(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "impute_method": 1,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "impute_method must be a string" in str(excinfo.value)

    def test_sia_linkage(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "sia_linkage": True,
        }
        cfg = config.Config(new)
        assert cfg.get("sia_linkage") is True

    def test_sia_linkage_not_bool(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "sia_linkage": 1,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "sia_linkage must be a boolean" in str(excinfo.value)

    def test_formula_file(self, input_file, output_file, clean_dir):
        formula_file = clean_dir / "formula.txt"
        formula_file.touch()
        formula_file = str(formula_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "formula_file": formula_file,
        }
        cfg = config.Config(new)
        assert cfg.get("formula_file") == formula_file

    def test_formula_file_none(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "formula_file": None,
        }
        cfg = config.Config(new)
        assert cfg.get("formula_file") is None

    def test_formula_file_not_str(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "formula_file": 1,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file must be a string" in str(excinfo.value)

    def test_formula_file_not_txt(self, input_file, output_file, clean_dir):
        formula_file = clean_dir / "formula.csv"
        formula_file.touch()
        formula_file = str(formula_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "formula_file": formula_file,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file must be a TXT file" in str(excinfo.value)

    def test_formula_file_not_exist(self, input_file, output_file, clean_dir):
        formula_file = clean_dir / "formula.txt"
        formula_file = str(formula_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "formula_file": formula_file,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file does not exist" in str(excinfo.value)

    def test_filter_invalid_traits(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "filter_invalid_traits": True,
        }
        cfg = config.Config(new)
        assert cfg.get("filter_invalid_traits") is True

    def test_group_file(self, input_file, output_file, clean_dir):
        group_file = clean_dir / "group.csv"
        group_file.touch()
        group_file = str(group_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "group_file": group_file,
        }
        cfg = config.Config(new)
        assert cfg.get("group_file") == group_file

    def test_group_file_not_str(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "group_file": 1,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file must be a string" in str(excinfo.value)

    def test_group_file_not_csv(self, input_file, output_file, clean_dir):
        group_file = clean_dir / "group.txt"
        group_file.touch()
        group_file = str(group_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "group_file": group_file,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file must be a CSV file" in str(excinfo.value)

    def test_group_file_not_exist(self, input_file, output_file, clean_dir):
        group_file = clean_dir / "group.csv"
        group_file = str(group_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "group_file": group_file,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file does not exist" in str(excinfo.value)

    def test_invalid_key(self, input_file, output_file):
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "invalid_key": 1,
        }
        with pytest.raises(KeyError) as excinfo:
            config.Config(new)
        assert "Invalid config key" in str(excinfo.value)
