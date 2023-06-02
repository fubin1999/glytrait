from typing import Optional

from glytrait.io import read_input, write_output
from glytrait.trait import load_formulas, build_meta_property_table, calcu_derived_trait, \
    calcu_direct_trait


def run_workflow(
    input_file: str,
    output_file: str,
    sia_linkage: bool = False,
    user_formula_file: Optional[str] = None,
) -> None:
    """Run the workflow."""
    glycans, abund_df = read_input(input_file)
    formulas = load_formulas(user_formula_file)
    if not sia_linkage:
        formulas = [f for f in formulas if f.sia_linkage is False]
    meta_prop_df = build_meta_property_table(abund_df.columns, glycans, sia_linkage)
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)
    direct_traits = calcu_direct_trait(abund_df)
    write_output(output_file, derived_traits, direct_traits, meta_prop_df, formulas)
