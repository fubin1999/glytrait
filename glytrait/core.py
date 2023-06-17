import pandas as pd

from glytrait.config import Config
from glytrait.io import read_input, write_output, read_group
from glytrait.preprocessing import preprocess_pipeline
from glytrait.stats import auto_hypothesis_test
from glytrait.trait import (
    load_formulas,
    build_meta_property_table,
    calcu_derived_trait,
    filter_derived_trait,
)


def run_workflow(config: Config) -> None:
    """Run the whole GlyTrait workflow.

    Args:
        config: Configuration for the workflow.
    """
    glycans, abund_df = read_input(config.get("input_file"))
    abund_df = preprocess_pipeline(
        abund_df, config.get("filter_glycan_max_na"), config.get("impute_method")
    )
    formulas = load_formulas(config.get("formula_file"))
    if not config.get("sia_linkage"):
        formulas = [f for f in formulas if f.sia_linkage is False]
    meta_prop_df = build_meta_property_table(
        abund_df.columns, glycans, config.get("sia_linkage")
    )
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)
    if config.get("filter_invalid_traits"):
        derived_traits = filter_derived_trait(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]
    if config.get("group_file") is not None:
        groups = read_group(config.get("group_file"))
        combined_traits = pd.concat([abund_df, derived_traits], axis=1)
        hypo_result = auto_hypothesis_test(combined_traits, groups)
    else:
        groups = None
        hypo_result = None
    write_output(
        config, derived_traits, abund_df, meta_prop_df, formulas, groups, hypo_result
    )
