from typing import Optional

import pandas as pd

from glytrait.io import read_input, write_output, read_group
from glytrait.stats import auto_hypothesis_test
from glytrait.trait import (
    load_formulas,
    build_meta_property_table,
    calcu_derived_trait,
    calcu_direct_trait,
    filter_derived_trait,
)


def run_workflow(
    input_file: str,
    output_file: str,
    sia_linkage: bool = False,
    user_formula_file: Optional[str] = None,
    filter_invalid: bool = True,
    group_file: Optional[str] = None,
) -> None:
    """Run the workflow."""
    glycans, abund_df = read_input(input_file)
    formulas = load_formulas(user_formula_file)
    if not sia_linkage:
        formulas = [f for f in formulas if f.sia_linkage is False]
    meta_prop_df = build_meta_property_table(abund_df.columns, glycans, sia_linkage)
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)
    if filter_invalid:
        derived_traits = filter_derived_trait(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]
    direct_traits = calcu_direct_trait(abund_df)
    if group_file is not None:
        groups = read_group(group_file)
        combined_traits = pd.concat([direct_traits, derived_traits], axis=1)
        hypo_result = auto_hypothesis_test(combined_traits, groups)
    else:
        hypo_result = None
    write_output(
        output_file, derived_traits, direct_traits, meta_prop_df, formulas, hypo_result
    )
