import collections

import pandas as pd
import pytest
from attrs import define

from glytrait import workflow as gw
from glytrait.config import default_config as glytrait_default_config, Config


@pytest.fixture
def default_config(monkeypatch) -> Config:
    monkeypatch.setattr(Config, "validators", [])
    new = dict(
        input_file="input_file",
        output_file="output_file",
        mode="structure",
        filter_glycan_max_na=0.5,
        impute_method="min",
        corr_threshold=1.0,
        corr_method="pearson",
        sia_linkage=False,
        formula_file=None,
        post_filtering=True,
        group_file=None,
        structure_file=None,
        database=None,
    )
    config_dict = glytrait_default_config.copy()
    config_dict.update(new)
    config = Config(config_dict)
    return config


class TestWorkflowState:
    @pytest.fixture
    def state(self):
        return gw.WorkflowState({"key": "value"})

    def test_get(self, state):
        assert state.get("key") == "value"

    def test_set(self, state):
        state.set("key", "new_value")
        assert state.get("key") == "new_value"

    def test_set_new(self, state):
        state.set("new_key", "new_value")
        assert state.get("new_key") == "new_value"

    def test_clear(self, state):
        state.clear()
        assert state.get("key") is None


class TestWorkflowStep:
    def test_skip(self, default_config):
        @define
        class TestStep(gw.WorkflowStep):
            def _check_needed(self):
                return False

            def _execute(self):
                return {"key": "value"}

        state = gw.WorkflowState()
        step = TestStep(default_config, state)
        step.run()
        assert state.get("key") is None

    def test_run(self, default_config):
        @define
        class TestStep(gw.WorkflowStep):
            def _execute(self):
                return {"key": "value"}

        state = gw.WorkflowState()
        step = TestStep(default_config, state)
        step.run()
        assert state.get("key") == "value"


@pytest.fixture
def state():
    return gw.WorkflowState()


@pytest.fixture
def input_df_basic():
    return pd.DataFrame(
        {
            "Composition": ["comp1", "comp2", "comp3"],
            "Structure": ["str1", "str2", "str3"],
            "Sample1": [1, 2, 3],
            "Sample2": [4, 5, 6],
            "Sample3": [7, 8, 9],
        }
    )


@pytest.fixture
def input_df_no_struc_col(input_df_basic):
    return input_df_basic.drop("Structure", axis=1)


class TestReadInputFileStep:
    def test_basic(self, mocker, default_config, state, input_df_basic):
        mocker.patch(
            "glytrait.workflow.pd.read_csv", return_value=input_df_basic, autospec=True
        )
        step = gw.ReadInputFileStep(default_config, state)
        step.run()
        assert state.get("input_df").equals(input_df_basic)


class TestCheckInputFilesStep:
    def test_basic(self, mocker, default_config, state, input_df_basic):
        state.set("input_df", input_df_basic)
        check_input_mock = mocker.patch(
            "glytrait.workflow.check_input_file", autospec=True
        )
        step = gw.CheckInputFileStep(default_config, state)
        step.run()
        check_input_mock.assert_called_once_with(input_df_basic)
        assert state.get("has_struc_col") is True

    def test_no_struc_col(self, mocker, default_config, state, input_df_no_struc_col):
        state.set("input_df", input_df_no_struc_col)
        check_input_mock = mocker.patch(
            "glytrait.workflow.check_input_file", autospec=True
        )
        step = gw.CheckInputFileStep(default_config, state)
        step.run()
        check_input_mock.assert_called_once_with(input_df_no_struc_col)
        assert state.get("has_struc_col") is False


class TestCheckHasStructureStep:
    def test_has_struc_col(self, default_config, state):
        state.set("has_struc_col", True)
        step = gw.CheckHasStructureStep(default_config, state)
        step.run()

    def test_has_db(self, default_config, state):
        state.set("has_struc_col", False)
        default_config.set("database", "db")
        step = gw.CheckHasStructureStep(default_config, state)
        step.run()

    def test_has_struc_file(self, default_config, state):
        state.set("has_struc_col", False)
        default_config.set("structure_file", "structure_file")
        step = gw.CheckHasStructureStep(default_config, state)
        step.run()

    def test_no_struc_info(self, default_config, state):
        state.set("has_struc_col", False)
        step = gw.CheckHasStructureStep(default_config, state)
        with pytest.raises(gw.InputError) as excinfo:
            step.run()
        msg = (
            "Must provide either structure_file or database name when the input file "
            "does not have a 'Structure' column."
        )
        assert msg in str(excinfo.value)


class TestLoadAbundanceTableStep:
    @pytest.mark.parametrize("has_struc_col", [True, False])
    def test_with_struc_col(
        self,
        default_config,
        state,
        has_struc_col,
        input_df_basic,
        input_df_no_struc_col,
    ):
        if has_struc_col:
            input_df = input_df_basic
        else:
            input_df = input_df_no_struc_col
        state.set("input_df", input_df)
        step = gw.LoadAbundanceTableStep(default_config, state)
        step.run()
        result_df = pd.DataFrame(
            {
                "comp1": [1, 4, 7],
                "comp2": [2, 5, 8],
                "comp3": [3, 6, 9],
            },
            index=["Sample1", "Sample2", "Sample3"],
        )
        assert state.get("abund_df").equals(result_df)


class TestLoadGlycansStep:
    @pytest.fixture
    def load_comp_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.load_compositions", return_value="glycans", autospec=True
        )

    @pytest.fixture
    def read_struc_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.read_structure_file",
            return_value="glycans",
            autospec=True,
        )

    @pytest.fixture
    def load_default_struc_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.load_default_structures",
            return_value="glycans",
            autospec=True,
        )

    @pytest.fixture
    def load_glycans_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.load_glycans", return_value="glycans", autospec=True
        )

    def test_from_struc_col(
        self,
        mocker,
        default_config,
        state,
        input_df_basic,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_basic)
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        assert state.get("glycans") == "glycans"
        load_glycans_mock.assert_called_once()
        load_comp_mock.assert_not_called()
        load_default_struc_mock.assert_not_called()
        read_struc_mock.assert_not_called()

    def test_from_db(
        self,
        mocker,
        default_config,
        state,
        input_df_no_struc_col,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_no_struc_col)
        default_config.set("database", "db")
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        assert state.get("glycans") == "glycans"
        load_glycans_mock.assert_not_called()
        load_comp_mock.assert_not_called()
        load_default_struc_mock.assert_called_once()
        read_struc_mock.assert_not_called()

    def test_from_struc_file(
        self,
        mocker,
        default_config,
        state,
        input_df_no_struc_col,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_no_struc_col)
        default_config.set("structure_file", "structure_file")
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        assert state.get("glycans") == "glycans"
        load_glycans_mock.assert_not_called()
        load_comp_mock.assert_not_called()
        load_default_struc_mock.assert_not_called()
        read_struc_mock.assert_called_once()

    def test_struc_col_and_db(
        self,
        mocker,
        default_config,
        state,
        input_df_basic,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_basic)
        default_config.set("database", "db")
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        load_glycans_mock.assert_called_once()
        load_comp_mock.assert_not_called()
        load_default_struc_mock.assert_not_called()
        read_struc_mock.assert_not_called()

    def test_struc_col_and_struc_file(
        self,
        mocker,
        default_config,
        state,
        input_df_basic,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_basic)
        default_config.set("structure_file", "structure_file")
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        load_glycans_mock.assert_called_once()
        load_comp_mock.assert_not_called()
        load_default_struc_mock.assert_not_called()
        read_struc_mock.assert_not_called()

    def test_comp(
        self,
        mocker,
        default_config,
        state,
        input_df_no_struc_col,
        load_comp_mock,
        load_default_struc_mock,
        read_struc_mock,
        load_glycans_mock,
    ):
        state.set("input_df", input_df_no_struc_col)
        default_config.set("mode", "composition")
        step = gw.LoadGlycansStep(default_config, state)
        step.run()
        assert state.get("glycans") == "glycans"
        load_glycans_mock.assert_not_called()
        load_comp_mock.assert_called_once()
        load_default_struc_mock.assert_not_called()
        read_struc_mock.assert_not_called()


class TestLoadGroupStep:
    @pytest.fixture
    def abund_df(self):
        return pd.DataFrame(
            {
                "comp1": [1, 2, 3, 4, 5, 6],
                "comp2": [7, 8, 9, 10, 11, 12],
                "comp3": [13, 14, 15, 16, 17, 18],
            },
            index=["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"],
        )

    @pytest.fixture
    def groups(self) -> pd.Series:
        samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"]
        return pd.Series(["A", "A", "A", "B", "B", "B"], name="group", index=samples)

    @pytest.fixture
    def read_group_mock(self, mocker, groups):
        return mocker.patch(
            "glytrait.workflow.read_group_file", return_value=groups, autospec=True
        )

    def test_basic(self, default_config, groups, read_group_mock, state, abund_df):
        default_config.set("group_file", "group_file")
        state.set("abund_df", abund_df)
        step = gw.LoadGroupStep(default_config, state)
        step.run()
        assert state.get("groups").equals(groups)
        read_group_mock.assert_called_once_with("group_file")

    def test_no_group_file(self, default_config, read_group_mock, state, abund_df):
        state.set("abund_df", abund_df)
        step = gw.LoadGroupStep(default_config, state)
        step.run()
        assert state.get("groups") is None
        read_group_mock.assert_not_called()

    def test_diff_samples(self, default_config, read_group_mock, state, abund_df):
        default_config.set("group_file", "group_file")
        state.set("abund_df", abund_df.drop("Sample6", axis=0))
        step = gw.LoadGroupStep(default_config, state)
        with pytest.raises(gw.InputError) as excinfo:
            step.run()
        msg = "The group file must have the same samples as the input file."
        assert msg in str(excinfo.value)
        read_group_mock.assert_called_once_with("group_file")

    def test_order(self, default_config, groups, read_group_mock, state, abund_df):
        new_index = ["Sample2", "Sample3", "Sample1", "Sample4", "Sample5", "Sample6"]
        groups_ordered = groups.reindex(new_index)
        read_group_mock.return_value = groups_ordered
        default_config.set("group_file", "group_file")
        state.set("abund_df", abund_df)
        step = gw.LoadGroupStep(default_config, state)
        step.run()
        assert state.get("groups").equals(groups)
        read_group_mock.assert_called_once_with("group_file")

    def test_only_1_group(
        self, default_config, groups, read_group_mock, state, abund_df
    ):
        groups = pd.Series(["A"] * 6, name="group", index=groups.index)
        read_group_mock.return_value = groups
        default_config.set("group_file", "group_file")
        state.set("abund_df", abund_df)
        step = gw.LoadGroupStep(default_config, state)
        with pytest.raises(gw.InputError) as excinfo:
            step.run()
        msg = "The group file must have at least 2 groups."
        assert msg in str(excinfo.value)
        read_group_mock.assert_called_once_with("group_file")

    def test_too_few_samples_in_groups(
        self, default_config, groups, read_group_mock, state, abund_df
    ):
        groups = groups.drop(["Sample1", "Sample4"])
        read_group_mock.return_value = groups
        abund_df = abund_df.drop(["Sample1", "Sample4"], axis=0)
        state.set("abund_df", abund_df)
        default_config.set("group_file", "group_file")
        step = gw.LoadGroupStep(default_config, state)
        with pytest.raises(gw.InputError) as excinfo:
            step.run()
        msg = "The following groups have less than 3 samples: A, B."
        assert msg in str(excinfo.value)
        read_group_mock.assert_called_once_with("group_file")


class TestLoadFormulasStep:
    @define
    class Formula:
        name: str
        sia_linkage: bool

    @pytest.fixture
    def formulas(self) -> list:
        return [
            self.Formula("formula1", False),
            self.Formula("formula2", True),
        ]

    @pytest.fixture
    def formulas_no_sia(self) -> list:
        return [self.Formula("formula1", False)]

    @pytest.fixture
    def load_formulas_mock(self, mocker, formulas):
        return mocker.patch(
            "glytrait.workflow.load_formulas", return_value=formulas, autospec=True
        )

    def test_basic(
        self,
        default_config,
        state,
        formulas,
        formulas_no_sia,
        load_formulas_mock,
    ):
        default_config.set("formula_file", "formula_file")
        step = gw.LoadFormulasStep(default_config, state)
        step.run()
        assert state.get("formulas") == formulas_no_sia
        load_formulas_mock.assert_called_once_with("structure", "formula_file")

    def test_sia_linkage(self, default_config, state, formulas, load_formulas_mock):
        default_config.set("formula_file", "formula_file")
        default_config.set("sia_linkage", True)
        step = gw.LoadFormulasStep(default_config, state)
        step.run()
        assert state.get("formulas") == formulas
        load_formulas_mock.assert_called_once_with("structure", "formula_file")

    def test_comp_mode(
        self, default_config, state, formulas, formulas_no_sia, load_formulas_mock
    ):
        default_config.set("formula_file", "formula_file")
        default_config.set("mode", "composition")
        step = gw.LoadFormulasStep(default_config, state)
        step.run()
        assert state.get("formulas") == formulas_no_sia
        load_formulas_mock.assert_called_once_with("composition", "formula_file")


class TestPreprocessStep:
    def test_basic(self, mocker, default_config, state):
        Glycan = collections.namedtuple("Glycan", ["name"])
        glycans = [Glycan("glycan1"), Glycan("glycan2")]

        state.set("abund_df", "abund_df")
        state.set("glycans", glycans)

        abund_df_processed = mocker.Mock()
        abund_df_processed.columns = ["glycan1"]
        preprocess_mock = mocker.patch(
            "glytrait.workflow.preprocess_pipeline",
            return_value=abund_df_processed,
            autospec=True,
        )

        step = gw.PreprocessStep(default_config, state)
        step.run()

        assert state.get("abund_df") == abund_df_processed
        assert state.get("glycans") == [Glycan("glycan1")]
        preprocess_mock.assert_called_once_with("abund_df", 0.5, "min")


class TestCalcTraitStep:
    @pytest.fixture
    def calcu_meta_prop_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.build_meta_property_table",
            return_value="meta_prop",
            autospec=True,
        )

    @pytest.fixture
    def calcu_traits_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.calcu_derived_trait",
            return_value="derived_traits",
            autospec=True,
        )

    @pytest.fixture
    def abund_df(self):
        return pd.DataFrame(
            {
                "comp1": [1, 2, 3, 4, 5, 6],
                "comp2": [7, 8, 9, 10, 11, 12],
                "comp3": [13, 14, 15, 16, 17, 18],
            },
            index=["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"],
        )

    def test_basic(
        self, default_config, state, abund_df, calcu_meta_prop_mock, calcu_traits_mock
    ):
        state.set("abund_df", abund_df)
        state.set("glycans", "glycans")
        state.set("formulas", "formulas")

        step = gw.CalcTraitStep(default_config, state)
        step.run()

        assert state.get("meta_property_df") == "meta_prop"
        assert state.get("derived_trait_df") == "derived_traits"

        calcu_meta_prop_mock.assert_called_once_with(
            ["comp1", "comp2", "comp3"], "glycans", "structure", False
        )
        calcu_traits_mock.assert_called_once_with(abund_df, "meta_prop", "formulas")


class TestPostFilteringStep:
    @pytest.fixture
    def abund_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "comp1": [1, 2, 3],
                "comp2": [4, 5, 6],
                "comp3": [7, 8, 9],
            },
            index=["Sample1", "Sample2", "Sample3"],
        )

    @pytest.fixture
    def filter_invalid_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.filter_invalid",
            return_value=("formulas_1", "traits_1"),
            autospec=True,
        )

    @pytest.fixture
    def filter_collinear_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.filter_colinearity",
            return_value=("formulas_2", "traits_2"),
            autospec=True,
        )

    @pytest.fixture
    def state_ok(self, state, abund_df):
        state.set("traits_filtered", False)
        state.set("abund_df", abund_df)
        state.set("formulas", "formulas_0")
        state.set("derived_trait_df", "traits_0")
        return state

    def test_basic(
        self,
        default_config,
        state_ok,
        abund_df,
        filter_invalid_mock,
        filter_collinear_mock,
    ):
        step = gw.PostFilteringStep(default_config, state_ok)
        step.run()

        assert state_ok.get("formulas") == "formulas_2"
        assert state_ok.get("derived_trait_df") == "traits_2"
        assert state_ok.get("traits_filtered") is True
        filter_invalid_mock.assert_called_once_with("formulas_0", "traits_0")
        filter_collinear_mock.assert_called_once_with(
            "formulas_1", "traits_1", 1.0, "pearson"
        )

    def test_skip(
        self, default_config, state_ok, filter_invalid_mock, filter_collinear_mock
    ):
        default_config.set("post_filtering", False)

        step = gw.PostFilteringStep(default_config, state_ok)
        step.run()

        assert state_ok.get("formulas") == "formulas_0"
        assert state_ok.get("derived_trait_df") == "traits_0"
        assert state_ok.get("traits_filtered") is False
        filter_invalid_mock.assert_not_called()
        filter_collinear_mock.assert_not_called()

    def test_too_few_samples(
        self,
        default_config,
        state_ok,
        filter_invalid_mock,
        filter_collinear_mock,
        abund_df,
    ):
        abund_df = abund_df.drop("Sample3", axis=0)
        state_ok.set("abund_df", abund_df)

        step = gw.PostFilteringStep(default_config, state_ok)
        step.run()

        assert state_ok.get("formulas") == "formulas_0"
        assert state_ok.get("derived_trait_df") == "traits_0"
        assert state_ok.get("traits_filtered") is False
        filter_invalid_mock.assert_not_called()
        filter_collinear_mock.assert_not_called()


class TestAnalysisStep:
    @pytest.fixture
    def abund_df(self):
        return pd.DataFrame(
            {
                "comp1": [1] * 9,
                "comp2": [1] * 9,
                "comp3": [1] * 9,
            },
            index=[f"Sample{i}" for i in range(1, 10)],
        )

    @pytest.fixture
    def traits(self):
        return pd.DataFrame(
            {
                "trait1": [1] * 9,
                "trait2": [1] * 9,
                "trait3": [1] * 9,
            },
            index=[f"Sample{i}" for i in range(1, 10)],
        )

    @pytest.fixture
    def combined_df(self, abund_df, traits):
        return pd.concat([abund_df, traits], axis=1)

    @pytest.fixture
    def groups(self) -> pd.Series:
        samples = [f"Sample{i}" for i in range(1, 10)]
        return pd.Series(["A"] * 4 + ["B"] * 5, name="group", index=samples)

    @pytest.fixture
    def differential_analysis_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.differential_analysis",
            return_value="hypo_test_result",
            autospec=True,
        )

    @pytest.fixture
    def roc_mock(self, mocker):
        return mocker.patch(
            "glytrait.workflow.calcu_roc_auc",
            return_value="roc_result",
            autospec=True,
        )

    @pytest.fixture
    def state_ok(self, state, abund_df, traits, groups):
        state.set("abund_df", abund_df)
        state.set("derived_trait_df", traits)
        state.set("groups", groups)
        state.set("traits_filtered", True)
        return state

    def test_basic(self, default_config, state_ok, differential_analysis_mock, roc_mock):
        step = gw.AnalysisStep(default_config, state_ok)
        step.run()

        assert state_ok.get("univariate_result") == "hypo_test_result"
        assert state_ok.get("roc_result") == "roc_result"

        differential_analysis_mock.assert_called_once()
        roc_mock.assert_called_once()

    def test_no_groups(self, default_config, state_ok, differential_analysis_mock, roc_mock):
        state_ok.set("groups", None)
        step = gw.AnalysisStep(default_config, state_ok)
        step.run()

        assert state_ok.get("univariate_result") is None
        assert state_ok.get("roc_result") is None

        differential_analysis_mock.assert_not_called()
        roc_mock.assert_not_called()

    def test_not_filtered(self, default_config, state_ok, differential_analysis_mock, roc_mock):
        state_ok.set("traits_filtered", False)
        step = gw.AnalysisStep(default_config, state_ok)
        step.run()

        assert state_ok.get("univariate_result") is None
        assert state_ok.get("roc_result") is None

        differential_analysis_mock.assert_not_called()
        roc_mock.assert_not_called()

    def test_3_groups(self, default_config, state_ok, differential_analysis_mock, roc_mock):
        groups = pd.Series(
            ["A"] * 3 + ["B"] * 3 + ["C"] * 3,
            name="group",
            index=state_ok.get("groups").index,
        )
        state_ok.set("groups", groups)

        step = gw.AnalysisStep(default_config, state_ok)
        step.run()

        assert state_ok.get("univariate_result") == "hypo_test_result"
        assert state_ok.get("roc_result") is None

        differential_analysis_mock.assert_called_once()
        roc_mock.assert_not_called()


class TestWriteOutputStep:
    @pytest.fixture
    def state_ok(self, state):
        state.set("abund_df", "abund_df")
        state.set("glycans", "glycans")
        state.set("formulas", "formulas")
        state.set("groups", "groups")
        state.set("meta_property_df", "meta_property_df")
        state.set("derived_trait_df", "derived_trait_df")
        state.set("traits_filtered", True)
        state.set("univariate_result", "univariate_result")
        state.set("roc_result", "roc_result")
        return state

    def test_basic(self, mocker, default_config, state_ok):
        write_output_mock = mocker.patch(
            "glytrait.workflow.write_output",
            autospec=True,
        )
        step = gw.WriteOutputStep(default_config, state_ok)
        step.run()

        write_output_mock.assert_called_once_with(
            default_config,
            "derived_trait_df",
            "abund_df",
            "meta_property_df",
            "formulas",
            "groups",
            "univariate_result",
            "roc_result",
        )
