from collections import namedtuple
from typing import Any, ClassVar

import pandas as pd
import pytest
from attrs import define, field

import glytrait.api
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


@pytest.mark.skip("This module will be refactored.")
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
        ]
        glytrait_api.export_all.assert_called_once_with(to_export, output_path)

        glytrait_api.load_formulas_from_file.assert_called_once_with(
            "structure", None, False
        )
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

        glytrait_api.load_formulas_from_file.assert_called_once_with(
            "structure", "custom_formula_file", False
        )
        glytrait_api.calcu_derived_trait.assert_called_once_with(
            "abundance_table", "meta_property_table", "formulas"
        )


@define
class SimpleWorkflow(glytrait_api._Workflow):
    _all_steps: ClassVar = ["step1", "step2", "step3"]
    _data_dict: ClassVar = {
        "step1": ["data1"],
        "step2": ["data2"],
        "step3": ["data3_1", "data3_2"],
    }

    @glytrait_api._step
    def step1(self, return_value=None):
        if return_value:
            return {"data1": return_value}
        return {"data1": "abcde"}

    @glytrait_api._step
    def step2(self):
        return {"data2": 12345}

    @glytrait_api._step
    def step3(self):
        return {"data3_1": "xyz", "data3_2": 67890}


class TestWorkflow:

    def test_call_step1(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        assert workflow.get_data("data1") == "abcde"

    def test_call_step1_and_step2(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        workflow.step2()
        assert workflow.get_data("data1") == "abcde"
        assert workflow.get_data("data2") == 12345

    def test_call_step1_step2_and_step3(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        workflow.step2()
        workflow.step3()
        assert workflow.get_data("data1") == "abcde"
        assert workflow.get_data("data2") == 12345
        assert workflow.get_data("data3_1") == "xyz"
        assert workflow.get_data("data3_2") == 67890

    def test_recall_step1_no_force(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        with pytest.raises(glytrait_api.InvalidOperationOrderError):
            workflow.step1()

    def test_recall_step1_with_force(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        workflow.step1(return_value="new_value", force=True)
        assert workflow.get_data("data1") == "new_value"

    def test_reset_data(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        workflow.step2()
        workflow.step1(force=True)
        assert workflow.get_data("data1") == "abcde"
        with pytest.raises(glytrait_api.MissingDataError):
            workflow.get_data("data2")

    def test_call_step2_without_step1(self):
        workflow = SimpleWorkflow()
        with pytest.raises(glytrait_api.InvalidOperationOrderError):
            workflow.step2()

    def test_call_step3_without_step2(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        with pytest.raises(glytrait_api.InvalidOperationOrderError):
            workflow.step3()

    def test_call_step3_without_step1(self):
        workflow = SimpleWorkflow()
        with pytest.raises(glytrait_api.InvalidOperationOrderError):
            workflow.step3()

    def test_current_step(self):
        workflow = SimpleWorkflow()
        assert workflow.current_step == "__START__"

        workflow.step1()
        assert workflow.current_step == "step1"

        workflow.step2()
        assert workflow.current_step == "step2"

        workflow.step3()
        assert workflow.current_step == "step3"

        workflow.step1(force=True)
        assert workflow.current_step == "step1"

    def test_reset(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        workflow.step2()
        workflow.reset()
        assert workflow.current_step == "__START__"
        with pytest.raises(glytrait_api.MissingDataError):
            workflow.get_data("data1")
        with pytest.raises(glytrait_api.MissingDataError):
            workflow.get_data("data2")

    def test_get_data_not_generated_yet(self):
        workflow = SimpleWorkflow()
        with pytest.raises(glytrait_api.MissingDataError):
            workflow.get_data("data1")

    def test_get_data_not_exist(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        with pytest.raises(KeyError):
            workflow.get_data("not_exist")


@define
class FakeInputData:
    abundance_table: pd.DataFrame
    glycans: dict[str, str]
    groups: pd.Series


class TestExperiment:

    @pytest.fixture
    def abundance_table(self):
        return pd.DataFrame(
            {
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )

    @pytest.fixture
    def glycans(self):
        return {"G1": "Glycan1", "G2": "Glycan2", "G3": "Glycan3"}

    @pytest.fixture
    def groups(self):
        return pd.Series(
            data=["A", "B", "A"],
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
            name="Group",
        )

    @pytest.fixture
    def input_data(self, abundance_table, glycans, groups):
        return FakeInputData(
            abundance_table=abundance_table, glycans=glycans, groups=groups
        )

    @pytest.fixture
    def exp(self, mocker, input_data):
        mocker.patch("glytrait.api.Experiment.reset")
        return glytrait_api.Experiment(input_data)

    def test_abundance_table_getter(self, exp, abundance_table):
        assert exp.abundance_table.equals(abundance_table)

    def test_abundance_table_setter(self, exp):
        exp.abundance_table = pd.DataFrame()
        assert exp.abundance_table.empty
        exp.reset.assert_called_once()

    def test_glycans_getter(self, exp, glycans):
        assert exp.glycans == glycans

    def test_glycans_setter(self, exp):
        exp.glycans = {}
        assert exp.glycans == {}
        exp.reset.assert_called_once()

    def test_groups_getter(self, exp, groups):
        assert exp.groups.equals(groups)

    def test_groups_setter(self, exp):
        exp.groups = pd.Series()
        assert exp.groups.empty
        exp.reset.assert_called_once()

    @pytest.mark.parametrize("filter", [0.0, 1.0])
    @pytest.mark.parametrize("impute_method", ["zero", "min"])
    def test_preprocess(
        self, mocker, input_data, abundance_table, filter, impute_method
    ):
        mocker.patch("glytrait.api.preprocess", return_value="result")
        exp = glytrait_api.Experiment(input_data)
        exp.preprocess(filter, impute_method)
        assert exp.processed_abundance_table == "result"
        called_args = glytrait_api.preprocess.call_args.kwargs
        assert called_args["data"].equals(abundance_table)
        assert called_args["filter_max_na"] == filter
        assert called_args["impute_method"] == impute_method

    def test_extract_meta_properties(self, mocker, input_data, abundance_table):
        mocker.patch("glytrait.api.build_meta_property_table", return_value="result")
        input_data.glycans["G4"] = "Glycan4"
        exp = glytrait_api.Experiment(input_data)
        exp._data["processed_abundance_table"] = abundance_table
        exp._current_step = "preprocess"
        exp.extract_meta_properties()
        assert exp.meta_property_table == "result"
        glytrait_api.build_meta_property_table.assert_called_once_with(
            {"G1": "Glycan1", "G2": "Glycan2", "G3": "Glycan3"}, "structure", False
        )

    @pytest.fixture
    def patch_for_derive_traits(self, mocker):
        mocker.patch("glytrait.api.load_default_formulas", return_value="formulas")
        mocker.patch("glytrait.api.parse_formulas")
        mocker.patch("glytrait.api.calcu_derived_trait", return_value="trait_table")

    @pytest.fixture
    def exp_for_derive_traits(self, exp):
        exp._data["processed_abundance_table"] = "abund_table"
        exp._data["meta_property_table"] = "mp_table"
        exp._current_step = "extract_meta_properties"
        return exp

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_default_formulas(self, exp_for_derive_traits):
        exp_for_derive_traits.derive_traits()
        assert exp_for_derive_traits.derived_trait_table == "trait_table"
        assert exp_for_derive_traits.get_data("formulas") == "formulas"
        glytrait.api.parse_formulas.assert_not_called()
        glytrait.api.calcu_derived_trait.assert_called_once()

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_provided_formulas(self, mocker, exp_for_derive_traits):
        FakeFormula = namedtuple("FakeFormula", "name sia_linkage")
        fake_formulas = [FakeFormula("F1", False), FakeFormula("F2", False)]
        mocker.patch("glytrait.api.parse_formulas", return_value=fake_formulas)
        exp_for_derive_traits.derive_traits(["expr1", "expr2"])
        assert exp_for_derive_traits.derived_trait_table == "trait_table"
        assert exp_for_derive_traits.get_data("formulas") == fake_formulas
        glytrait.api.load_default_formulas.assert_not_called()
        glytrait.api.calcu_derived_trait.assert_called_once()

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_provided_sia_formulas_no_sia_mode(
        self, mocker, exp_for_derive_traits
    ):
        FakeFormula = namedtuple("FakeFormula", "name sia_linkage")
        fake_formulas = [FakeFormula("F1", True), FakeFormula("F2", False)]
        mocker.patch("glytrait.api.parse_formulas", return_value=fake_formulas)
        with pytest.raises(ValueError):
            exp_for_derive_traits.derive_traits(["expr1", "expr2"])

    @pytest.mark.parametrize("corr_threshold", [1.0, 0.5])
    def test_post_filter(self, mocker, exp, corr_threshold):
        mocker.patch("glytrait.api.post_filter", return_value="filtered")
        exp._data["formulas"] = "formulas"
        exp._data["derived_trait_table"] = "trait_table"
        exp._current_step = "derive_traits"
        exp.post_filter(corr_threshold)
        glytrait.api.post_filter.assert_called_once_with(
            formulas="formulas",
            trait_df="trait_table",
            threshold=corr_threshold,
            method="pearson",
        )
        assert exp.filtered_derived_trait_table == "filtered"

    @pytest.fixture
    def patch_for_diff_analysis(self, mocker):
        mocker.patch("glytrait.api.t_test", return_value="t_test_result")
        mocker.patch("glytrait.api.anova", return_value="anova_result")

    @pytest.fixture
    def exp_for_diff_analysis(self, exp):
        exp._data["filtered_derived_trait_table"] = "trait_table"
        exp._current_step = "post_filter"
        return exp

    @pytest.mark.usefixtures("patch_for_diff_analysis")
    def test_diff_analysis(self, exp_for_diff_analysis):
        exp_for_diff_analysis.input_data.groups = None
        with pytest.raises(glytrait.api.MissingDataError):
            exp_for_diff_analysis.diff_analysis()
