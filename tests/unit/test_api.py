from typing import Any

import pytest
from attrs import define, field

from glytrait.api import _Config, GlyTrait
from glytrait import api as glytrait_api


class TestConfig:
    @pytest.fixture
    def config_dict(self) -> dict:
        return dict(
            mode="structure",
            filter_max_na=0.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=0.9,
            sia_linkage=False,
            custom_formula_file=None,
        )

    def test_init(self, config_dict):
        config = _Config(**config_dict)
        assert config.mode == "structure"
        assert config.filter_max_na == 0.0
        assert config.impute_method == "zero"
        assert config.post_filtering is True
        assert config.correlation_threshold == 0.9
        assert config.sia_linkage is False
        assert config.custom_formula_file is None

    def test_default_config(self):
        config = _Config()
        assert config.mode == "structure"
        assert config.filter_max_na == 1.0
        assert config.impute_method == "zero"
        assert config.post_filtering is True
        assert config.correlation_threshold == 1.0
        assert config.sia_linkage is False
        assert config.custom_formula_file is None

    @pytest.mark.parametrize("mode", ["invalid", 1, None])
    def test_invalid_mode(self, config_dict, mode):
        config_dict["mode"] = mode
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("filter_max_na", ["invalid", -1, 2, None])
    def test_invalid_filter_max_na(self, config_dict, filter_max_na):
        config_dict["filter_max_na"] = filter_max_na
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("impute_method", ["invalid", 1, None])
    def test_invalid_impute_method(self, config_dict, impute_method):
        config_dict["impute_method"] = impute_method
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("post_filtering", ["invalid", 1, None])
    def test_invalid_post_filtering(self, config_dict, post_filtering):
        config_dict["post_filtering"] = post_filtering
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("correlation_threshold", ["invalid", -2, -0.5, 2, None])
    def test_invalid_correlation_threshold(self, config_dict, correlation_threshold):
        config_dict["correlation_threshold"] = correlation_threshold
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("sia_linkage", ["invalid", 1, None])
    def test_invalid_sia_linkage(self, config_dict, sia_linkage):
        config_dict["sia_linkage"] = sia_linkage
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("custom_formula_file", [1])
    def test_invalid_custom_formula_file(self, config_dict, custom_formula_file):
        config_dict["custom_formula_file"] = custom_formula_file
        with pytest.raises(ValueError):
            _Config(**config_dict)


@define
class MockInputData:
    """Mock input data."""

    abundance_table = field()
    glycans = field()
    groups = field()


class TestGlyTrait:
    @pytest.fixture(autouse=True)
    def patch_export_all(self, mocker):
        mocker.patch("glytrait.api.export_all", return_value=None)

    @pytest.fixture(autouse=True)
    def patch_load_formulas(self, mocker):
        mocker.patch("glytrait.api.load_formulas", return_value="formulas")

    @pytest.fixture
    def mock_input_data(self):
        return MockInputData(
            abundance_table="abundance_table",
            glycans="glycans",
            groups="groups",
        )

    @pytest.fixture(autouse=True)
    def patch_load_input_data(self, mocker, mock_input_data):
        mocker.patch(
            "glytrait.api.load_input_data_from_csv", return_value=mock_input_data
        )

    @pytest.fixture(autouse=True)
    def patch_build_meta_property_table(self, mocker):
        mocker.patch(
            "glytrait.api.build_meta_property_table", return_value="meta_property_table"
        )

    @pytest.fixture(autouse=True)
    def patch_post_filter(self, mocker):
        mocker.patch(
            "glytrait.api.post_filter", return_value="filtered_derived_trait_table"
        )

    @pytest.fixture(autouse=True)
    def patch_calcu_preprocess(self, mocker):
        mocker.patch("glytrait.api.preprocess", return_value=None)

    @pytest.fixture(autouse=True)
    def patch_calcu_derived_trait(self, mocker):
        mocker.patch(
            "glytrait.api.calcu_derived_trait", return_value="derived_trait_table"
        )

    @pytest.fixture(autouse=True)
    def patch_differential_analysis(self, mocker):
        mocker.patch("glytrait.api.differential_analysis", return_value="diff_result")

    @pytest.fixture
    def default_config_dict(self) -> dict[str, Any]:
        return dict(
            mode="structure",
            filter_max_na=0.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=0.9,
            sia_linkage=False,
            custom_formula_file=None,
        )

    def test_init(self, default_config_dict):
        glytrait = GlyTrait(**default_config_dict)
        assert glytrait._config.mode == "structure"
        assert glytrait._config.filter_max_na == 0.0
        assert glytrait._config.impute_method == "zero"
        assert glytrait._config.post_filtering is True
        assert glytrait._config.correlation_threshold == 0.9
        assert glytrait._config.sia_linkage is False
        assert glytrait._config.custom_formula_file is None

    def test_run(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        to_export = [
            ("meta_properties.csv", "meta_property_table"),
            ("derived_traits.csv", "derived_trait_table"),
            ("glycan_abundance_processed.csv", "abundance_table"),
            ("derived_traits_filtered.csv", "filtered_derived_trait_table"),
            ("differential_analysis.csv", "diff_result"),
        ]
        glytrait_api.export_all.assert_called_once_with(to_export, output_path)

        glytrait_api.load_formulas.assert_called_once_with("structure", None, False)
        glytrait_api.load_input_data_from_csv.assert_called_once_with(
            abundance_file="abundance_file",
            glycan_file="glycan_file",
            group_file="group_file",
            mode="structure",
        )
        glytrait_api.preprocess.assert_called_once_with(mock_input_data, 0.0, "zero")
        glytrait_api.build_meta_property_table.assert_called_once_with(
            "glycans", "structure", False
        )
        glytrait_api.calcu_derived_trait.assert_called_once_with(
            "abundance_table", "meta_property_table", "formulas"
        )
        glytrait_api.post_filter.assert_called_once_with(
            formulas="formulas",
            trait_df="derived_trait_table",
            threshold=0.9,
            method="pearson",
        )
        glytrait_api.differential_analysis.assert_called_once_with(
            "filtered_derived_trait_table", "groups"
        )

    def test_run_mode(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        default_config_dict["mode"] = "composition"
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.load_input_data_from_csv.assert_called_once_with(
            abundance_file="abundance_file",
            glycan_file="glycan_file",
            group_file="group_file",
            mode="composition",
        )
        glytrait_api.build_meta_property_table.assert_called_once_with(
            "glycans", "composition", False
        )

    def test_run_filter_max_na(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        default_config_dict["filter_max_na"] = 0.5
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.preprocess.assert_called_once_with(mock_input_data, 0.5, "zero")

    def test_run_impute_method(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        default_config_dict["impute_method"] = "mean"
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.preprocess.assert_called_once_with(mock_input_data, 0.0, "mean")

    def test_run_post_filtering(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        default_config_dict["post_filtering"] = False
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.post_filter.assert_not_called()
        glytrait_api.export_all.assert_called_once_with(
            [
                ("meta_properties.csv", "meta_property_table"),
                ("derived_traits.csv", "derived_trait_table"),
                ("glycan_abundance_processed.csv", "abundance_table"),
            ],
            output_path,
        )

    def test_run_correlation_threshold(
        self, default_config_dict, mock_input_data, clean_dir
    ):
        output_path = str(clean_dir)
        default_config_dict["correlation_threshold"] = 0.5
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.post_filter.assert_called_once_with(
            formulas="formulas",
            trait_df="derived_trait_table",
            threshold=0.5,
            method="pearson",
        )

    def test_run_sia_linkage(self, default_config_dict, mock_input_data, clean_dir):
        output_path = str(clean_dir)
        default_config_dict["sia_linkage"] = True
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.build_meta_property_table.assert_called_once_with(
            "glycans", "structure", True
        )

    def test_run_custom_formula_file(
        self, default_config_dict, mock_input_data, clean_dir
    ):
        output_path = str(clean_dir)
        default_config_dict["custom_formula_file"] = "custom_formula_file"
        glytrait = GlyTrait(**default_config_dict)
        glytrait.run(
            output_path,
            "abundance_file",
            "glycan_file",
            "group_file",
        )

        glytrait_api.load_formulas.assert_called_once_with(
            "structure", "custom_formula_file", False
        )
        glytrait_api.calcu_derived_trait.assert_called_once_with(
            "abundance_table", "meta_property_table", "formulas"
        )
