import pytest

import glytrait.workflow
from glytrait.config import default_config as glytrait_default_config, Config


@pytest.fixture
def default_config(monkeypatch) -> Config:
    monkeypatch.setattr(glytrait.workflow.Config, "validators", [])
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


class TestWorkflow:

    def test_run_basic(self, mocker, default_config):
        # ----- Mock functions -----
        auto_hypothesis_test_mock = mocker.create_autospec("glytrait.workflow.auto_hypothesis_test")
        calcu_roc_auc_mock = mocker.create_autospec("glytrait.workflow.calcu_roc_auc")
        load_formulas_mock = mocker.create_autospec("glytrait.workflow.load_formulas")
        load_glycans_mock = mocker.create_autospec("glytrait.workflow.load_glycans")
        load_compositions_mock = mocker.create_autospec("glytrait.workflow.load_compositions")
        read_input_mock = mocker.create_autospec("glytrait.workflow.read_input")
        write_output_mock = mocker.create_autospec("glytrait.workflow.write_output")
        read_group_mock = mocker.create_autospec("glytrait.workflow.read_group")
        read_structure_mock = mocker.create_autospec("glytrait.workflow.read_structure")
        load_default_structures_mock = mocker.create_autospec("glytrait.workflow.load_default_structures")
        build_meta_property_table_mock = mocker.create_autospec("glytrait.workflow.build_meta_property_table")
        preprocess_pipeline_mock = mocker.create_autospec("glytrait.workflow.preprocess_pipeline")
        calcu_derived_trait_mock = mocker.create_autospec("glytrait.workflow.calcu_derived_trait")
        filter_invalid_mock = mocker.create_autospec("glytrait.workflow.filter_invalid")
        filter_colinearity_mock = mocker.create_autospec("glytrait.workflow.filter_colinearity")
        # ----- Mock functions END -----

        # ----- Set return values -----
        # auto_hypothesis_test_mock is not called
        # calcu_roc_auc_mock is not called

        formulas_mock = mocker.Mock()
        load_formulas_mock.return_value = formulas_mock

        glycans_mock = mocker.Mock()
        load_glycans_mock.return_value = glycans_mock

        # load_compositions_mock is not called

        comps_mock = mocker.Mock()
        strucs_mock = mocker.Mock()
        abund_table_mock = mocker.Mock()
        read_input_mock.return_value = comps_mock, strucs_mock, abund_table_mock

        # write_output_mock has no return value
        # read_group_mock is not called
        # read_structure_mock is not called
        # load_default_structures_mock is not called

        meta_prop_table_mock = mocker.Mock()
        build_meta_property_table_mock.return_value = meta_prop_table_mock

        preprocessed_mock = mocker.Mock()
        preprocess_pipeline_mock.return_value = preprocessed_mock

        derived_trait_mock = mocker.Mock()
        calcu_derived_trait_mock.return_value = derived_trait_mock

        filtered_formulas_mock1 = mocker.Mock()
        filtered_trait_table_mock1 = mocker.Mock()
        filter_invalid_mock.return_value = filtered_formulas_mock1, filtered_trait_table_mock1

        filtered_formulas_mock2 = mocker.Mock()
        filtered_trait_table_mock2 = mocker.Mock()
        filter_colinearity_mock.return_value = filtered_formulas_mock2, filtered_trait_table_mock2
        # ----- Set return values END -----
        