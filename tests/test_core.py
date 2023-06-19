from unittest.mock import Mock

import pytest

import glytrait.core


@pytest.mark.parametrize(
    "sia_linkage, filter_invalid, has_structure, database, structure_file, group_file",
    # sia_linkage means what the sia_linkage option is set
    # filter_invalid means what the filter_invalid_traits option is set
    # has_structure means whether the input file has a structure column
    # database means whether the input file is a database
    # structure_file means whether the structure file is provided
    # group_file means whether the group file is provided

    # Note that has_structure, database and structure_file are mutually exclusive
    [
        (True, True, True, None, None, None),
        (False, True, True, None, None, None),
        (True, False, True, None, None, None),
        (False, True, False, None, "structure_file", None),
        (False, True, False, "database", None, None),
        (False, True, True, None, None, "group_file"),
    ],
)
def test_run_workflow(
    mocker,
    sia_linkage,
    filter_invalid,
    has_structure,
    database,
    structure_file,
    group_file,
):
    raw_abund_df_mock = Mock(name="raw_abund_df_mock")
    raw_abund_df_mock.columns = ["col1", "col2"]
    abund_df_mock = Mock(name="abund_df_mock")
    abund_df_mock.columns = ["col1", "col2"]
    meta_prop_df_mock = Mock(name="meta_prop_df_mock")
    derived_trait_df_mock = Mock(name="derived_trait_df_mock")
    derived_trait_df_filtered_mock = Mock(name="derived_trait_df_filtered_mock")
    derived_trait_df_filtered_mock.columns = ["trait1"]
    group_series_mock = Mock(name="group_series_mock") if group_file else None
    hypo_test_result_mock = Mock(name="hypo_test_result_mock") if group_file else None
    combined_traits_mock = Mock(name="combined_traits_mock")

    formula_1_mock = Mock()
    formula_1_mock.sia_linkage = False
    formula_1_mock.name = "trait1"
    formula_2_mock = Mock()
    formula_2_mock.sia_linkage = True
    formula_2_mock.name = "trait2"
    formula_3_mock = Mock()
    formula_3_mock.sia_linkage = False
    formula_3_mock.name = "trait3"
    formulas_mock = [formula_1_mock, formula_2_mock, formula_3_mock]

    if has_structure:
        read_input_return_value = ("glycans", raw_abund_df_mock)
    else:
        read_input_return_value = (None, raw_abund_df_mock)
    read_input_mock = mocker.patch(
        "glytrait.core.read_input", return_value=read_input_return_value
    )
    read_structure_mock = mocker.patch(
        "glytrait.core.read_structure", return_value="glycans"
    )
    load_default_mock = mocker.patch(
        "glytrait.core.load_default", return_value="glycans"
    )
    preprocess_pipeline_mock = mocker.patch(
        "glytrait.core.preprocess_pipeline", return_value=abund_df_mock
    )
    load_formulas_mock = mocker.patch(
        "glytrait.core.load_formulas", return_value=formulas_mock
    )
    build_meta_property_table_mock = mocker.patch(
        "glytrait.core.build_meta_property_table", return_value=meta_prop_df_mock
    )
    calcu_derived_traits_mock = mocker.patch(
        "glytrait.core.calcu_derived_trait", return_value=derived_trait_df_mock
    )
    write_output_mock = mocker.patch("glytrait.core.write_output")
    filter_derived_trait_mock = mocker.patch(
        "glytrait.core.filter_derived_trait",
        return_value=derived_trait_df_filtered_mock,
    )
    read_group_mock = mocker.patch(
        "glytrait.core.read_group", return_value=group_series_mock
    )
    hypo_test_mock = mocker.patch(
        "glytrait.core.auto_hypothesis_test", return_value=hypo_test_result_mock
    )
    pd_concat_mock = mocker.patch("pandas.concat", return_value=combined_traits_mock)

    config = dict(
        input_file="input_file",
        output_file="output_file",
        filter_glycan_max_na=0.5,
        impute_method="min",
        sia_linkage=sia_linkage,
        formula_file="user_formula_file",
        filter_invalid_traits=filter_invalid,
        group_file=group_file,
        structure_file=structure_file,
        database=database,
    )
    glytrait.core.run_workflow(config)

    read_input_mock.assert_called_once_with("input_file")
    if has_structure:
        read_structure_mock.assert_not_called()
        load_default_mock.assert_not_called()
    elif structure_file:
        read_structure_mock.assert_called_once_with(
            "structure_file", abund_df_mock.columns
        )
        load_default_mock.assert_not_called()
    else:  # database
        read_structure_mock.assert_not_called()
        load_default_mock.assert_called_once_with("database", abund_df_mock.columns)
    preprocess_pipeline_mock.assert_called_once_with(raw_abund_df_mock, 0.5, "min")
    load_formulas_mock.assert_called_once_with("user_formula_file")
    build_meta_property_table_mock.assert_called_once_with(
        ["col1", "col2"], "glycans", sia_linkage
    )
    if filter_invalid:
        filter_derived_trait_mock.assert_called_once_with(derived_trait_df_mock)
    else:
        filter_derived_trait_mock.assert_not_called()

    if sia_linkage:
        calcu_derived_traits_mock.assert_called_once_with(
            abund_df_mock, meta_prop_df_mock, formulas_mock
        )
    else:
        calcu_derived_traits_mock.assert_called_once_with(
            abund_df_mock, meta_prop_df_mock, [formula_1_mock, formula_3_mock]
        )

    if filter_invalid:
        formulas_mock = [formula_1_mock]
    elif sia_linkage:
        formulas_mock = [formula_1_mock, formula_2_mock, formula_3_mock]
    else:
        formulas_mock = [formula_1_mock, formula_3_mock]

    if filter_invalid:
        derived_trait_df = derived_trait_df_filtered_mock
    else:
        derived_trait_df = derived_trait_df_mock

    if group_file:
        read_group_mock.assert_called_once_with("group_file")
        hypo_test_mock.assert_called_once_with(combined_traits_mock, group_series_mock)
    else:
        read_group_mock.assert_not_called()
        hypo_test_mock.assert_not_called()

    write_output_mock.assert_called_once_with(
        config,
        derived_trait_df,
        abund_df_mock,
        meta_prop_df_mock,
        formulas_mock,
        group_series_mock,
        hypo_test_result_mock,
    )
