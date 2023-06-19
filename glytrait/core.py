import pandas as pd

from glytrait.config import Config
from glytrait.io import read_input, write_output, read_group, read_structure
from glytrait.preprocessing import preprocess_pipeline
from glytrait.stats import auto_hypothesis_test
from glytrait.trait import (
    load_formulas,
    build_meta_property_table,
    calcu_derived_trait,
    filter_derived_trait,
)
from glytrait.database import load_default


def run_workflow(config: Config) -> None:
    """Run the whole GlyTrait workflow.

    Args:
        config: Configuration for the workflow.
    """
    # Load glycans and the abundance table
    glycans, abund_df = read_input(config.get("input_file"))
    if glycans is None:
        config.update({"_has_struc_col": False})
        if (struc_file := config.get("structure_file")) is not None:
            glycans = read_structure(struc_file, abund_df.columns)
        else:  # config.get("database") is not None
            glycans = load_default(config.get("database"), abund_df.columns)
    else:
        config.update({"_has_struc_col": True})

    # Preprocess the abundance table
    abund_df = preprocess_pipeline(
        abund_df, config.get("filter_glycan_max_na"), config.get("impute_method")
    )

    # Load formulas
    formulas = load_formulas(config.get("formula_file"))
    if not config.get("sia_linkage"):
        formulas = [f for f in formulas if f.sia_linkage is False]

    # Calculate derived traits
    meta_prop_df = build_meta_property_table(
        abund_df.columns, glycans, config.get("sia_linkage")
    )
    derived_traits = calcu_derived_trait(abund_df, meta_prop_df, formulas)

    # Filter derived traits
    if config.get("filter_invalid_traits"):
        derived_traits = filter_derived_trait(derived_traits)
        formulas = [f for f in formulas if f.name in derived_traits.columns]

    # Perform hypothesis test if group file is provided
    if config.get("group_file") is not None:
        groups = read_group(config.get("group_file"))
        combined_traits = pd.concat([abund_df, derived_traits], axis=1)
        hypo_result = auto_hypothesis_test(combined_traits, groups)
    else:
        groups = None
        hypo_result = None

    # Write output
    write_output(
        config, derived_traits, abund_df, meta_prop_df, formulas, groups, hypo_result
    )
