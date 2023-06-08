from unittest.mock import Mock

import pytest

import glytrait.core


@pytest.mark.parametrize(
    "sia_linkage, filter_invalid",
    [
        (True, True),
        (False, True),
        (True, False),
    ],
)
def test_run_workflow(mocker, sia_linkage, filter_invalid):
    abund_df_mock = Mock()
    abund_df_mock.columns = ["col1", "col2"]
    meta_prop_df_mock = Mock()
    derived_trait_df_mock = Mock()
    derived_trait_df_filtered_mock = Mock()
    derived_trait_df_filtered_mock.columns = ["trait1"]
    direct_trait_df_mock = Mock()

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

    read_input_mock = mocker.patch(
        "glytrait.core.read_input", return_value=("glycans", abund_df_mock)
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
    calcu_direct_traits_mock = mocker.patch(
        "glytrait.core.calcu_direct_trait", return_value=direct_trait_df_mock
    )
    write_output_mock = mocker.patch("glytrait.core.write_output")
    filter_derived_trait_mock = mocker.patch(
        "glytrait.core.filter_derived_trait",
        return_value=derived_trait_df_filtered_mock,
    )

    glytrait.core.run_workflow(
        "input_file", "output_file", sia_linkage, "user_formula_file", filter_invalid
    )

    read_input_mock.assert_called_once_with("input_file")
    load_formulas_mock.assert_called_once_with("user_formula_file")
    build_meta_property_table_mock.assert_called_once_with(
        ["col1", "col2"], "glycans", sia_linkage
    )
    calcu_direct_traits_mock.assert_called_once_with(abund_df_mock)
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
    write_output_mock.assert_called_once_with(
        "output_file",
        derived_trait_df,
        direct_trait_df_mock,
        meta_prop_df_mock,
        formulas_mock,
    )
