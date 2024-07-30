from collections import namedtuple
from typing import ClassVar

import pandas as pd
import pytest
from attrs import define

from glytrait import api


@define
class SimpleWorkflow(api._Workflow):
    _all_steps: ClassVar = ["step1", "step2", "step3"]
    _data_dict: ClassVar = {
        "step1": ["data1"],
        "step2": ["data2"],
        "step3": ["data3_1", "data3_2"],
    }

    @api._step
    def step1(self, return_value=None):
        if return_value:
            return {"data1": return_value}
        return {"data1": "abcde"}

    @api._step
    def step2(self):
        return {"data2": 12345}

    @api._step
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
        with pytest.raises(api.InvalidOperationOrderError):
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
        with pytest.raises(api.MissingDataError):
            workflow.get_data("data2")

    def test_call_step2_without_step1(self):
        workflow = SimpleWorkflow()
        with pytest.raises(api.InvalidOperationOrderError):
            workflow.step2()

    def test_call_step3_without_step2(self):
        workflow = SimpleWorkflow()
        workflow.step1()
        with pytest.raises(api.InvalidOperationOrderError):
            workflow.step3()

    def test_call_step3_without_step1(self):
        workflow = SimpleWorkflow()
        with pytest.raises(api.InvalidOperationOrderError):
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
        with pytest.raises(api.MissingDataError):
            workflow.get_data("data1")
        with pytest.raises(api.MissingDataError):
            workflow.get_data("data2")

    def test_get_data_not_generated_yet(self):
        workflow = SimpleWorkflow()
        with pytest.raises(api.MissingDataError):
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

    # The `Experiment` api needs files to be passed in as arguments.
    # So, the following fixtures store the data in files and return the file paths.
    @pytest.fixture
    def abundance_file(self, abundance_table, tmp_path):
        file_path = tmp_path / "abundance.csv"
        abundance_table.to_csv(file_path)
        return str(file_path)

    @pytest.fixture
    def glycans_file(self, glycans, tmp_path):
        file_path = tmp_path / "glycans.csv"
        pd.DataFrame(glycans.items(), columns=["GlycanID", "Structure"]).to_csv(file_path)
        return str(file_path)

    @pytest.fixture
    def groups_file(self, groups, tmp_path):
        file_path = tmp_path / "groups.csv"
        groups.reset_index().to_csv(file_path)
        return str(file_path)

    @pytest.fixture
    def exp(self, mocker, abundance_file, groups_file, glycans_file, glycans):
        mocker.patch("glytrait.api.Experiment.reset")
        mocker.patch("glytrait.api.load_glycans", return_value=glycans)
        exp = api.Experiment(
            abundance_file=abundance_file,
            glycan_file=glycans_file,
            group_file=groups_file,
        )
        return exp

    def test_abundance_table_getter(self, exp, abundance_table):
        assert exp.abundance_table.equals(abundance_table)

    def test_glycans_getter(self, exp, glycans):
        assert exp.glycans == glycans

    def test_groups_getter(self, exp, groups):
        assert exp.groups.equals(groups)

    @pytest.mark.parametrize("filter", [0.0, 1.0])
    @pytest.mark.parametrize("impute_method", ["zero", "min"])
    def test_preprocess(
        self, mocker, monkeypatch, abundance_table, filter, impute_method, exp
    ):
        mocker.patch("glytrait.api.build_meta_property_table", return_value="mp_table")
        mocker.patch("glytrait.api.preprocess", return_value=abundance_table + 1)
        exp.preprocess(filter, impute_method)
        assert exp.processed_abundance_table.equals(abundance_table + 1)
        assert exp.meta_property_table == "mp_table"
        called_args = api.preprocess.call_args.kwargs
        assert called_args["data"].equals(abundance_table)
        assert called_args["filter_max_na"] == filter
        assert called_args["impute_method"] == impute_method

    @pytest.fixture
    def patch_for_derive_traits(self, mocker):
        mocker.patch("glytrait.api.load_default_formulas", return_value="formulas")
        mocker.patch("glytrait.api.calcu_derived_trait", return_value="trait_table")

    @pytest.fixture
    def exp_for_derive_traits(self, exp):
        exp._data["processed_abundance_table"] = "abund_table"
        exp._data["meta_property_table"] = "mp_table"
        exp._current_step = "preprocess"
        return exp

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_default_formulas(self, exp_for_derive_traits):
        exp_for_derive_traits.derive_traits()
        assert exp_for_derive_traits.derived_trait_table == "trait_table"
        assert exp_for_derive_traits.get_data("formulas") == "formulas"
        api.calcu_derived_trait.assert_called_once()

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_provided_formulas(self, exp_for_derive_traits):
        FakeFormula = namedtuple("FakeFormula", "name sia_linkage")
        fake_formulas = [FakeFormula("F1", False), FakeFormula("F2", False)]
        exp_for_derive_traits.derive_traits(fake_formulas)
        assert exp_for_derive_traits.derived_trait_table == "trait_table"
        assert exp_for_derive_traits.get_data("formulas") == fake_formulas
        api.load_default_formulas.assert_not_called()
        api.calcu_derived_trait.assert_called_once()

    @pytest.mark.usefixtures("patch_for_derive_traits")
    def test_derive_traits_with_provided_sia_formulas_no_sia_mode(
        self, exp_for_derive_traits
    ):
        FakeFormula = namedtuple("FakeFormula", "name sia_linkage")
        fake_formulas = [FakeFormula("F1", True), FakeFormula("F2", False)]
        with pytest.raises(ValueError):
            exp_for_derive_traits.derive_traits(fake_formulas)

    @pytest.mark.parametrize("corr_threshold", [1.0, 0.5])
    def test_post_filter(self, mocker, exp, corr_threshold):
        mocker.patch("glytrait.api.post_filter", return_value="filtered")
        exp._data["formulas"] = "formulas"
        exp._data["derived_trait_table"] = "trait_table"
        exp._current_step = "derive_traits"
        exp.post_filter(corr_threshold)
        api.post_filter.assert_called_once_with(
            formulas="formulas",
            trait_df="trait_table",
            threshold=corr_threshold,
            method="pearson",
        )
        assert exp.filtered_derived_trait_table == "filtered"

    @pytest.fixture
    def patch_for_diff_analysis(self, mocker):
        mocker.patch("glytrait.api.auto_test", return_value="t_test_result")

    @pytest.fixture
    def exp_for_diff_analysis(self, exp):
        exp._data["filtered_derived_trait_table"] = "trait_table"
        exp._current_step = "post_filter"
        return exp

    @pytest.mark.usefixtures("patch_for_diff_analysis")
    def test_diff_analysis_no_group(self, exp_for_diff_analysis):
        exp_for_diff_analysis.groups = None
        with pytest.raises(api.MissingDataError):
            exp_for_diff_analysis.diff_analysis()

    @pytest.mark.usefixtures("patch_for_diff_analysis")
    def test_diff_analysis(self, exp_for_diff_analysis):
        exp_for_diff_analysis.diff_analysis()
        assert exp_for_diff_analysis.diff_results == "t_test_result"
        api.auto_test.assert_called_once()

    @pytest.fixture
    def patch_for_run_workflow(self, mocker):
        mocker.patch("glytrait.api.Experiment.preprocess")
        mocker.patch("glytrait.api.Experiment.derive_traits")
        mocker.patch("glytrait.api.Experiment.post_filter")
        mocker.patch("glytrait.api.Experiment.diff_analysis")

    @pytest.mark.usefixtures("patch_for_run_workflow")
    def test_run_workflow_with_groups_no_args(self, exp):
        exp.run_workflow()

        exp.preprocess.assert_called_once_with(1.0, "zero")
        exp.derive_traits.assert_called_once_with(None)
        exp.post_filter.assert_called_once_with(1.0)
        exp.diff_analysis.assert_called_once()

    @pytest.mark.usefixtures("patch_for_run_workflow")
    def test_run_workflow_with_args(self, exp):
        exp.run_workflow(
            filter_max_na=0.5,
            impute_method="min",
            corr_threshold=0.5,
            formulas=["F1", "F2"],
        )

        exp.preprocess.assert_called_once_with(0.5, "min")
        exp.derive_traits.assert_called_once_with(["F1", "F2"])
        exp.post_filter.assert_called_once_with(0.5)
        exp.diff_analysis.assert_called_once()

    @pytest.mark.usefixtures("patch_for_run_workflow")
    def test_run_workflow_without_groups(self, exp):
        exp.groups = None
        exp.run_workflow()

        exp.preprocess.assert_called_once_with(1.0, "zero")
        exp.derive_traits.assert_called_once_with(None)
        exp.post_filter.assert_called_once_with(1.0)
        exp.diff_analysis.assert_not_called()

    def test_try_formulas_one(self, mocker, exp_for_derive_traits):
        trait_table = pd.DataFrame(
            {"F1": [1, 2, 3]},
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        mocker.patch("glytrait.api.calcu_derived_trait", return_value=trait_table)
        mocker.patch("glytrait.api.parse_formulas", return_value=["F1"])
        result = exp_for_derive_traits.try_formulas("F1")
        expected = pd.Series(
            [1, 2, 3], index=pd.Index(["S1", "S2", "S3"], name="Sample"), name="F1"
        )
        pd.testing.assert_series_equal(result, expected)
        api.calcu_derived_trait.assert_called_once()
        api.parse_formulas.assert_called_once_with(["F1"])

    def test_try_formulas_one_no_squeeze(self, mocker, exp_for_derive_traits):
        trait_table = pd.DataFrame(
            {"F1": [1, 2, 3]},
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        mocker.patch("glytrait.api.calcu_derived_trait", return_value=trait_table)
        mocker.patch("glytrait.api.parse_formulas", return_value=["F1"])
        result = exp_for_derive_traits.try_formulas("F1", squeeze=False)
        pd.testing.assert_frame_equal(result, trait_table)
        api.calcu_derived_trait.assert_called_once()
        api.parse_formulas.assert_called_once_with(["F1"])

    def test_try_formulas_multiple(self, mocker, exp_for_derive_traits):
        trait_table = pd.DataFrame(
            {
                "F1": [1, 2, 3],
                "F2": [4, 5, 6],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        mocker.patch("glytrait.api.calcu_derived_trait", return_value=trait_table)
        mocker.patch("glytrait.api.parse_formulas", return_value=["F1", "F2"])
        result = exp_for_derive_traits.try_formulas(["F1", "F2"])
        pd.testing.assert_frame_equal(result, trait_table)
        api.calcu_derived_trait.assert_called_once()
        api.parse_formulas.assert_called_once_with(["F1", "F2"])

    def test_try_formulas_no_preprocess(self, exp):
        with pytest.raises(api.InvalidOperationOrderError):
            exp.try_formulas("F1")
