from typing import Any

import pytest

import glytrait.workflow
from glytrait.config import default_config as glytrait_default_config


@pytest.fixture
def default_config(monkeypatch) -> dict[str, Any]:
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
        formula_file="user_formula_file",
        post_filtering=True,
        group_file=None,
        structure_file=None,
        database=None,
    )
    config = glytrait_default_config.copy()
    config.update(new)
    return config


@pytest.fixture
def formula_mocks(mocker) -> list:
    formula_1_mock = mocker.Mock()
    formula_1_mock.sia_linkage = False
    formula_1_mock.name = "trait1"
    formula_2_mock = mocker.Mock()
    formula_2_mock.sia_linkage = True
    formula_2_mock.name = "trait2"
    formula_3_mock = mocker.Mock()
    formula_3_mock.sia_linkage = False
    formula_3_mock.name = "trait3"
    return [formula_1_mock, formula_2_mock, formula_3_mock]


@pytest.mark.parametrize(
    "has_structure, mode, structure_file, database",
    [
        (True, "structure", None, None),
        (False, "structure", "structure_file", None),
        (False, "structure", None, "database"),
        (False, "composition", None, None),
    ],
)
def test_load_and_preprocess_data(
    mocker, default_config, mode, has_structure, structure_file, database
):
    config = default_config.copy()
    new = {
        "mode": mode,
        "structure_file": structure_file,
        "database": database,
    }
    config.update(new)

    comps = ["comp1", "comp2", "comp3"]
    strucs = ["struc1", "struc2", "struc3"]
    abund_df = mocker.Mock(name="abund_df_mock")
    abund_df.columns = comps
    glycans = ["glycan1", "glycan2", "glycan3"]
    strucs = None if not has_structure else strucs
    read_input_mock = mocker.patch(
        "glytrait.workflow.read_input",
        return_value=(comps, strucs, abund_df),
        autospec=True,
    )
    load_glycans_mock = mocker.patch(
        "glytrait.workflow.load_glycans",
        return_value=glycans,
        autospec=True,
    )
    load_default_mock = mocker.patch(
        "glytrait.workflow.load_default_structures",
        return_value=glycans,
        autospec=True,
    )
    read_structure_mock = mocker.patch(
        "glytrait.workflow.read_structure",
        return_value=glycans,
        autospec=True,
    )
    load_compositions_mock = mocker.patch(
        "glytrait.workflow.load_compositions",
        return_value=glycans,
        autospec=True,
    )
    preprocess_mock = mocker.patch(
        "glytrait.workflow.preprocess_pipeline",
        return_value=abund_df,
        autospec=True,
    )

    result = glytrait.workflow._load_and_preprocess_data(config)

    read_input_mock.assert_called_once_with("input_file")
    match (has_structure, mode, structure_file, database):
        case (True, "structure", None, None):
            load_glycans_mock.assert_called_once_with(strucs)
            load_default_mock.assert_not_called()
            read_structure_mock.assert_not_called()
            load_compositions_mock.assert_not_called()
        case (False, "structure", "structure_file", None):
            load_glycans_mock.assert_not_called()
            load_default_mock.assert_not_called()
            read_structure_mock.assert_called_once_with(structure_file, comps)
            load_compositions_mock.assert_not_called()
        case (False, "structure", None, "database"):
            load_glycans_mock.assert_not_called()
            load_default_mock.assert_called_once_with(database, comps)
            read_structure_mock.assert_not_called()
            load_compositions_mock.assert_not_called()
        case (False, "composition", None, None):
            load_glycans_mock.assert_not_called()
            load_default_mock.assert_not_called()
            read_structure_mock.assert_not_called()
            load_compositions_mock.assert_called_once_with(comps, sia_linkage=False)
    preprocess_mock.assert_called_once_with(abund_df, 0.5, "min")
    assert result == (glycans, abund_df)


def test_preprocess(mocker, default_config):
    config = default_config.copy()
    new = {
        "filter_glycan_max_na": 0.5,
        "impute_method": "min",
    }
    config.update(new)

    before_abund_df = mocker.Mock(name="before_abund_df_mock")
    before_abund_df.columns = ["comp1", "comp2", "comp3"]
    before_comps = ["comp1", "comp2", "comp3"]
    before_strucs = ["struc1", "struc2", "struc3"]
    after_abund_df = mocker.Mock(name="after_abund_df_mock")
    after_abund_df.columns = ["comp1", "comp3"]
    after_comps = ["comp1", "comp3"]
    after_strucs = ["struc1", "struc3"]
    mocker.patch(
        "glytrait.workflow.preprocess_pipeline",
        return_value=after_abund_df,
        autospec=True,
    )

    result_comps, result_strucs, result_abund_df = glytrait.workflow._preprocess(
        config, before_comps, before_strucs, before_abund_df
    )
    assert result_comps == after_comps
    assert result_strucs == after_strucs
    assert result_abund_df == after_abund_df


@pytest.mark.parametrize(
    "sia_linkage",
    [True, False],
)
def test_load_formulas(mocker, default_config, formula_mocks, sia_linkage):
    config = default_config.copy()
    new = {
        "sia_linkage": sia_linkage,
        "formula_file": "formula_file",
    }
    config.update(new)

    load_formulas_mock = mocker.patch(
        "glytrait.workflow.load_formulas",
        return_value=formula_mocks,
        autospec=True,
    )
    result = glytrait.workflow._load_formulas(config)
    load_formulas_mock.assert_called_once_with("structure", "formula_file")
    if sia_linkage:
        expected = formula_mocks
    else:
        expected = [formula_mocks[0], formula_mocks[2]]
    assert result == expected


@pytest.mark.parametrize(
    "sia_linkage",
    [True, False],
)
def test_calcu_derived_traits(mocker, default_config, sia_linkage):
    config = default_config.copy()
    config.update({"sia_linkage": sia_linkage})

    abund_df_mock = mocker.Mock(name="abund_df_mock")
    abund_df_mock.columns = mocker.Mock(name="columns_mock")
    build_meta_property_table = mocker.patch(
        "glytrait.workflow.build_meta_property_table",
        return_value="meta_prop_df",
        autospec=True,
    )
    calcu_derived_trait_mock = mocker.patch(
        "glytrait.workflow.calcu_derived_trait",
        return_value="derived_trait_df",
        autospec=True,
    )

    result = glytrait.workflow._calcu_derived_traits(
        config, abund_df_mock, "glycans", "formulas"
    )

    build_meta_property_table.assert_called_once_with(
        abund_df_mock.columns, "glycans", "structure", sia_linkage
    )
    calcu_derived_trait_mock.assert_called_once_with(
        abund_df_mock, "meta_prop_df", "formulas"
    )
    assert result == ("meta_prop_df", "derived_trait_df")


@pytest.mark.parametrize(
    "post_filtering, corr_threshold, corr_method",
    [
        (True, 0.5, "pearson"),
        (False, 0.5, "pearson"),
        (True, 0.5, "spearman"),
        (True, 1, "pearson"),
    ],
)
def test_filter_traits(mocker, default_config, post_filtering, corr_threshold, corr_method):
    config = default_config.copy()
    new = {
        "corr_threshold": corr_threshold,
        "corr_method": corr_method,
        "post_filtering": post_filtering,
    }
    config.update(new)

    filter_invalid_mock = mocker.patch(
        "glytrait.workflow.filter_invalid",
        return_value=("filtered_formulas_1", "filtered_trait_df_1"),
        autospec=True,
    )
    filter_colinearity_mock = mocker.patch(
        "glytrait.workflow.filter_colinearity",
        return_value=("filtered_formulas_2", "filtered_trait_df_2"),
        autospec=True,
    )

    result = glytrait.workflow._filter_traits(config, "trait_df", "formulas")
    if post_filtering:
        expected = ("filtered_formulas_2", "filtered_trait_df_2")
        filter_invalid_mock.assert_called_once_with("formulas", "trait_df")
        filter_colinearity_mock.assert_called_once_with(
            "filtered_formulas_1",
            "filtered_trait_df_1",
            corr_threshold,
            corr_method,
        )
    else:
        expected = ("formulas", "trait_df")
        filter_invalid_mock.assert_not_called()
        filter_colinearity_mock.assert_not_called()
    assert result == expected


@pytest.mark.parametrize(
    "has_groups, n_groups",
    [
        (True, 2),
        (True, 3),
        (False, 0),
    ],
)
def test_statistical_analysis(mocker, default_config, has_groups, n_groups):
    config = default_config.copy()
    if has_groups is True:
        config.update({"group_file": "group_file"})
    else:
        config.update({"group_file": None})

    group_mock = mocker.Mock()
    group_mock.nunique.return_value = n_groups
    get_groups_mock = mocker.patch(
        "glytrait.workflow._get_groups",
        return_value=group_mock if has_groups else None,
        autospec=True,
    )
    hypo_test_mock = mocker.patch(
        "glytrait.workflow.auto_hypothesis_test",
        return_value="hypo_test_result",
        autospec=True,
    )
    roc_mock = mocker.patch(
        "glytrait.workflow.calcu_roc_auc",
        return_value="roc_result",
        autospec=True,
    )

    if has_groups and n_groups == 2:
        return_value = (group_mock, "hypo_test_result", "roc_result")
    elif has_groups and n_groups > 2:
        return_value = (group_mock, "hypo_test_result", None)
    else:
        return_value = (None, None, None)
    assert glytrait.workflow._statistical_analysis(config, "trait_df") == return_value


def test_run_workflow(mocker):
    load_and_preprocess_data_mock = mocker.patch(
        "glytrait.workflow._load_and_preprocess_data",
        return_value=("glycans", "abund_df"),
        autospec=True,
    )
    load_formulas_mock = mocker.patch(
        "glytrait.workflow._load_formulas",
        return_value="formulas",
        autospec=True,
    )
    calcu_derived_traits_mock = mocker.patch(
        "glytrait.workflow._calcu_derived_traits",
        return_value=("meta_prop_df", "derived_trait_df"),
        autospec=True,
    )
    filter_traits_mock = mocker.patch(
        "glytrait.workflow._filter_traits",
        return_value=("formulas", "derived_trait_df"),
        autospec=True,
    )
    combine_data_mock = mocker.patch(
        "glytrait.workflow._combine_data",
        return_value="trait_df",
        autospec=True,
    )
    statistical_analysis_mock = mocker.patch(
        "glytrait.workflow._statistical_analysis",
        return_value=("groups", "hypo_test_result", "roc_result"),
        autospec=True,
    )
    write_output_mock = mocker.patch("glytrait.workflow.write_output", autospec=True)

    glytrait.workflow.run_workflow("config")

    load_and_preprocess_data_mock.assert_called_once_with("config")
    load_formulas_mock.assert_called_once_with("config")
    calcu_derived_traits_mock.assert_called_once_with(
        "config", "abund_df", "glycans", "formulas"
    )
    filter_traits_mock.assert_called_once_with(
        "config", "derived_trait_df", "formulas"
    )
    combine_data_mock.assert_called_once_with("abund_df", "derived_trait_df")
    statistical_analysis_mock.assert_called_once_with("config", "trait_df")
    write_output_mock.assert_called_once_with(
        "config",
        "derived_trait_df",
        "abund_df",
        "meta_prop_df",
        "formulas",
        "groups",
        "hypo_test_result",
        "roc_result",
    )
