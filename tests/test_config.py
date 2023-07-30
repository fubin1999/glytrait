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

    @pytest.fixture
    def base_config(self, input_file, output_file):
        return {
            "input_file": input_file,
            "output_file": output_file,
            "mode": "structure",
        }

    def test_init(self, base_config):
        cfg = config.Config(base_config)
        assert cfg.asdict() == config.default_config | base_config

    def test_init_without_input_file(self, output_file):
        new = {"output_file": output_file, "mode": "structure"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Missing config key: input_file" in str(excinfo.value)

    def test_input_file_not_csv(self, clean_dir, output_file):
        input_file = clean_dir / "test.txt"
        input_file.touch()
        input_file = str(input_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "mode": "structure",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a CSV file" in str(excinfo.value)

    def test_input_file_not_exists(self, clean_dir, output_file):
        input_file = str(clean_dir / "test.csv")
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "mode": "structure",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file does not exist" in str(excinfo.value)

    def test_input_not_a_file(self, clean_dir, output_file):
        input_file = clean_dir / "test.csv"
        input_file.mkdir()
        input_file = str(input_file)
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "mode": "structure",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a file, not directory" in str(excinfo.value)

    def test_input_file_not_str(self, output_file):
        new = {"input_file": 1, "output_file": output_file, "mode": "structure"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Input file must be a string" in str(excinfo.value)

    def test_output_file_not_xlsx(self, clean_dir, input_file):
        output_file = str(clean_dir / "test.txt")
        new = {
            "input_file": input_file,
            "output_file": output_file,
            "mode": "structure",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Output file must be a XLSX file" in str(excinfo.value)

    def test_filter_glycan_max_na(self, base_config):
        new = base_config | {"filter_glycan_max_na": 0.2}
        cfg = config.Config(new)
        assert cfg.get("filter_glycan_max_na") == 0.2

    def test_filter_glycan_max_na_out_of_range(self, base_config):
        new = base_config | {"filter_glycan_max_na": 1.2}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "filter_glycan_max_na must be between 0 and 1" in str(excinfo.value)

    def test_filter_glycan_max_na_str(self, base_config):
        new = base_config | {"filter_glycan_max_na": "0.2"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "filter_glycan_max_na must be a float" in str(excinfo.value)

    def test_filter_glycan_max_na_int(self, base_config):
        new = base_config | {"filter_glycan_max_na": 1}
        cfg = config.Config(new)
        assert cfg.get("filter_glycan_max_na") == 1

    def test_impute_method(self, base_config):
        new = base_config | {"impute_method": "min"}
        cfg = config.Config(new)
        assert cfg.get("impute_method") == "min"

    def test_impute_method_not_supported(self, base_config):
        new = base_config | {"impute_method": "foo"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "impute_method must be one of" in str(excinfo.value)

    def test_impute_method_not_str(self, base_config):
        new = base_config | {"impute_method": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "impute_method must be a string" in str(excinfo.value)

    def test_correlation_threshold(self, base_config):
        new = base_config | {"correlation_threshold": 0.5}
        cfg = config.Config(new)
        assert cfg.get("correlation_threshold") == 0.5

    def test_correlation_threshold_out_of_range(self, base_config):
        new = base_config | {"correlation_threshold": 1.2}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "correlation_threshold must be between 0 and 1" in str(excinfo.value)

    def test_correlation_threshold_str(self, base_config):
        new = base_config | {"correlation_threshold": "0.5"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "correlation_threshold must be a float" in str(excinfo.value)

    def test_sia_linkage(self, base_config):
        new = base_config | {"sia_linkage": True}
        cfg = config.Config(new)
        assert cfg.get("sia_linkage") is True

    def test_sia_linkage_not_bool(self, base_config):
        new = base_config | {"sia_linkage": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "sia_linkage must be a boolean" in str(excinfo.value)

    def test_formula_file(self, base_config, clean_dir):
        formula_file = clean_dir / "formula.txt"
        formula_file.touch()
        formula_file = str(formula_file)
        new = base_config | {"formula_file": formula_file}
        cfg = config.Config(new)
        assert cfg.get("formula_file") == formula_file

    def test_formula_file_none(self, base_config):
        new = base_config | {"formula_file": None}
        cfg = config.Config(new)
        assert cfg.get("formula_file") is None

    def test_formula_file_not_str(self, base_config):
        new = base_config | {"formula_file": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file must be a string" in str(excinfo.value)

    def test_formula_file_not_txt(self, base_config, clean_dir):
        formula_file = clean_dir / "formula.csv"
        formula_file.touch()
        formula_file = str(formula_file)
        new = base_config | {"formula_file": formula_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file must be a TXT file" in str(excinfo.value)

    def test_formula_file_not_exist(self, base_config, clean_dir):
        formula_file = clean_dir / "formula.txt"
        formula_file = str(formula_file)
        new = base_config | {"formula_file": formula_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Formula file does not exist" in str(excinfo.value)

    def test_filter_invalid_traits(self, base_config):
        new = base_config | {"filter_invalid_traits": True}
        cfg = config.Config(new)
        assert cfg.get("filter_invalid_traits") is True

    def test_filter_invalid_traits_not_bool(self, base_config):
        new = base_config | {"filter_invalid_traits": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "filter_invalid_traits must be a boolean" in str(excinfo.value)

    def test_group_file(self, base_config, clean_dir):
        group_file = clean_dir / "group.csv"
        group_file.touch()
        group_file = str(group_file)
        new = base_config | {"group_file": group_file}
        cfg = config.Config(new)
        assert cfg.get("group_file") == group_file

    def test_group_file_not_str(self, base_config):
        new = base_config | {"group_file": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file must be a string" in str(excinfo.value)

    def test_group_file_not_csv(self, base_config, clean_dir):
        group_file = clean_dir / "group.txt"
        group_file.touch()
        group_file = str(group_file)
        new = base_config | {"group_file": group_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file must be a CSV file" in str(excinfo.value)

    def test_group_file_not_exist(self, base_config, clean_dir):
        group_file = clean_dir / "group.csv"
        group_file = str(group_file)
        new = base_config | {"group_file": group_file}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Group file does not exist" in str(excinfo.value)

    def test_invalid_key(self, base_config):
        new = base_config | {"foo": "bar"}
        with pytest.raises(KeyError) as excinfo:
            config.Config(new)
        assert "Invalid config key" in str(excinfo.value)

    def test_database(self, base_config):
        new = base_config | {"database": "serum"}
        cfg = config.Config(new)
        assert cfg.get("database") == "serum"

    def test_database_unknown(self, base_config):
        new = base_config | {"database": "not_valid"}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Database must be one of: serum, IgG" in str(excinfo.value)

    def test_database_not_str(self, base_config):
        new = base_config | {"database": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Database must be a string" in str(excinfo.value)

    def test_update_two_invalid_values(self, base_config):
        """Check that when two values ara updated, the first one valid and the second one
        invalid, both are not updated."""
        cfg = config.Config(base_config)
        old_value = cfg.get("filter_glycan_max_na")
        new = {
            "filter_glycan_max_na": 0.2,  # valid
            "impute_method": "not_valid",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            cfg.update(new)
        assert "impute_method must be one of" in str(excinfo.value)
        assert cfg.get("filter_glycan_max_na") == old_value

    def test_both_database_and_structure_file(self, base_config, clean_dir):
        structure_file = clean_dir / "structure.csv"
        structure_file.touch()
        new = {
            "database": "serum",
            "structure_file": str(structure_file),
        }
        new = base_config | new
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Cannot provide both database and structure_file" in str(excinfo.value)

    def test_both_struc_col_and_database(self, base_config):
        new = base_config | {"database": "serum", "_has_struc_col": True}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        msg = "Cannot provide database when the input file already"
        assert msg in str(excinfo.value)

    def test_both_struc_col_and_structure_file(self, base_config, clean_dir):
        structure_file = clean_dir / "structure.csv"
        structure_file.touch()
        new = base_config | {
            "structure_file": str(structure_file),
            "_has_struc_col": True,
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        msg = "Cannot provide structure_file when the input file already"
        assert msg in str(excinfo.value)

    def test_has_struc_col_not_bool(self, base_config):
        new = base_config | {"_has_struc_col": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "_has_struc_col must be a boolean" in str(excinfo.value)

    def test_no_struc_col_struc_mode_no_database_no_structure_file(self, base_config):
        new = base_config | {"_has_struc_col": False}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Must provide either structure_file or database" in str(excinfo.value)

    @pytest.mark.parametrize("has_struc_col", [True, False])
    def test_comp_mode(self, base_config, has_struc_col):
        new = base_config | {"_has_struc_col": has_struc_col, "mode": "composition"}
        cfg = config.Config(new)
        assert cfg.get("_has_struc_col") is has_struc_col

    def test_comp_mode_database_provided(self, base_config):
        new = base_config | {
            "database": "serum",
            "mode": "composition",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Database is incompatible with 'Composition' mode." in str(excinfo.value)

    def test_comp_mode_structure_file_provided(self, base_config, clean_dir):
        structure_file = clean_dir / "structure.csv"
        structure_file.touch()
        new = base_config | {
            "structure_file": str(structure_file),
            "mode": "composition",
        }
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        msg = "Structure file is incompatible with 'Composition' mode."
        assert msg in str(excinfo.value)

    def test_structure_file_a_directory(self, base_config, clean_dir):
        new = base_config | {"structure_file": str(clean_dir)}
        assert clean_dir.is_dir()

    def test_structure_file_not_csv(self, base_config, clean_dir):
        structure_file = clean_dir / "structure.txt"
        structure_file.touch()
        new = base_config | {"structure_file": str(structure_file)}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Structure file must be a CSV file" in str(excinfo.value)

    def test_structure_file_not_exist_as_csv(self, base_config, clean_dir):
        structure_file = clean_dir / "structure.csv"
        new = base_config | {"structure_file": str(structure_file)}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Structure file does not exist" in str(excinfo.value)

    def test_structure_file_not_exist_as_dir(self, base_config, clean_dir):
        structure_file = clean_dir / "structure"
        new = base_config | {"structure_file": str(structure_file)}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Structure file does not exist" in str(excinfo.value)

    def test_structure_file_not_str(self, base_config):
        new = base_config | {"structure_file": 1}
        with pytest.raises(config.ConfigError) as excinfo:
            config.Config(new)
        assert "Structure file must be a string" in str(excinfo.value)
