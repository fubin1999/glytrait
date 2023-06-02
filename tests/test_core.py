from unittest.mock import Mock, MagicMock

import pytest

import glytrait.core


@pytest.mark.parametrize("sia_linkage", [True, False])
def test_run_workflow(mocker, sia_linkage):
    abund_df_mock = Mock()
    abund_df_mock.columns = ["col1", "col2"]
    formulas_mock = [MagicMock(), MagicMock(), MagicMock()]
    formulas_mock[0].sia_linkage = False
    formulas_mock[1].sia_linkage = True
    formulas_mock[2].sia_linkage = False

    read_input_mock = mocker.patch(
        "glytrait.core.read_input", return_value=("glycans", abund_df_mock)
    )
    load_formulas_mock = mocker.patch(
        "glytrait.core.load_formulas", return_value=formulas_mock
    )
    build_meta_property_table_mock = mocker.patch(
        "glytrait.core.build_meta_property_table", return_value="meta_prop_df"
    )
    calcu_derived_traits_mock = mocker.patch(
        "glytrait.core.calcu_derived_trait", return_value="derived_trait_df"
    )
    calcu_direct_traits_mock = mocker.patch(
        "glytrait.core.calcu_direct_trait", return_value="direct_trait_df"
    )
    write_output_mock = mocker.patch("glytrait.core.write_output")

    glytrait.core.run_workflow("input_file", "output_file", sia_linkage, "user_formula_file")

    read_input_mock.assert_called_once_with("input_file")
    load_formulas_mock.assert_called_once_with("user_formula_file")
    build_meta_property_table_mock.assert_called_once_with(
        ["col1", "col2"], "glycans", sia_linkage
    )

    if sia_linkage is False:
        formulas_mock.pop(1)
    calcu_derived_traits_mock.assert_called_once_with(
        abund_df_mock, "meta_prop_df", formulas_mock
    )
    write_output_mock.assert_called_once_with(
        "output_file",
        "derived_trait_df",
        "direct_trait_df",
        "meta_prop_df",
        formulas_mock,
    )
