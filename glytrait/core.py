from typing import Optional, Literal

import pandas as pd

from glytrait.io import read_input, write_output, read_group
from glytrait.preprocessing import preprocess_pipeline
from glytrait.stats import auto_hypothesis_test
from glytrait.trait import (
    load_formulas,
    build_meta_property_table,
    calcu_derived_trait,
    filter_derived_trait,
)


def run_workflow(
    input_file: str,
    output_file: str,
    *,
    filter_max_na: float = 0.5,
    impute_method: Literal["zero", "min", "lod", "mean", "median"] = "min",
    sia_linkage: bool = False,
    user_formula_file: Optional[str] = None,
    filter_invalid: bool = True,
    group_file: Optional[str] = None,
) -> None:
    """Run the whole GlyTrait workflow.

    Args:
        input_file (str): The input file path.
        output_file (str): The output file path.
        filter_max_na (float, optional): The maximum proportion of missing values allowed for a
            glycan. Defaults to 0.5.
        impute_method (Literal["zero", "min", "lod", "mean", "median"], optional): The imputation
            method. Can be "zero", "min", "1/5min", "mean", "median". Defaults to "min".
        sia_linkage (bool, optional): Whether to include sialic acid linkage information.
            Defaults to False.
        user_formula_file (Optional[str], optional): The user-defined formula file path. Defaults
            to None.
        filter_invalid (bool, optional): Whether to filter out invalid formulas. Defaults to True.
        group_file (Optional[str], optional): The group file path. Defaults to None.
    """
    glycans, abund_df = read_input(input_file)
    abund_df = preprocess_pipeline(abund_df, filter_max_na, impute_method)
    formulas = load_formulas(user_formula_file)
    if not sia_linkage:
        formulas = [f for f in formulas if f.sia_linkage is False]
    meta_prop_df = build_meta_property_table(abund_df.columns, glycans, sia_linkage)
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)
    if filter_invalid:
        derived_traits = filter_derived_trait(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]
    if group_file is not None:
        groups = read_group(group_file)
        combined_traits = pd.concat([abund_df, derived_traits], axis=1)
        hypo_result = auto_hypothesis_test(combined_traits, groups)
    else:
        hypo_result = None
    write_output(
        output_file, derived_traits, abund_df, meta_prop_df, formulas, hypo_result
    )
